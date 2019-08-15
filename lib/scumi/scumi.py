import sys
import os
import numpy as np
import scipy.io
import scipy.sparse
import numba
import random
import multiprocessing as mp
import subprocess
import cytoolz as toolz
import collections
from itertools import chain
import regex as re
import yaml
import logging
import time
import gzip
import pandas as pd
from functools import partial
from typing import NamedTuple
from pysam import AlignmentFile

from .util import compute_edit_distance, read_gene_map_from_gtf
from .fastq_io import read_fastq
from .barcode import ErrorBarcodeHash, ErrorBarcodeHashConstraint
from .estimate_cell_barcode import get_cell_whitelist

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s: %(levelname)s: %(message)s')
logger = logging.getLogger(__name__)


def format_fastq(*fastq, config, method, fastq_out, cb_count, 
                 num_thread=4, max_num_cell=1000000):
    """
    Merging fastq reads by putting the cell barcodes and UMI sequences
    to the headers of the cDNA reads

    :param config: the config file
    :param method: the library preparation protocol, e.g., can be one of
    10X, Drop-seq, InDrop, Seq-Well, CEL-seq2, sci-RNA-seq, SPLiT-seq,
    you can add protocol to the configure file easily
    by specifying the read structures.
    A template configuration file is provided in scumi/config.yaml
    :param fastq: input fastq files
    :param fastq_out: the output fastq file
    :param cb_count: an output file containing the # reads for each cell barcode
    :param num_thread: int
                     the number of cpu cores to use
    :param max_num_cell: int
                        the maximum number of cells
    """

    with open(config, 'r') as stream:
        config_dict = yaml.load(stream)

    config_dict = config_dict[method]
    num_read = config_dict['num_read']

    num_fastq = len(fastq)
    if num_fastq != num_read:
        logger.error(f'Error: the number of input fastq files {num_fastq} is different '
                     f'from the number of fastq files {num_read} detected in the config file')
        sys.exit(-1)

    read_regex_str, barcode_filter, read_regex_str_qual = \
        zip(*[_extract_input_read_template('read' + str(i), config_dict)
              for i in range(1, num_read + 1)])

    barcode_filter_dict = dict()
    for d in barcode_filter:
        barcode_filter_dict.update(d)

    read_template = _infer_read_template(read_regex_str)

    # select
    read_regex_list = [re.compile(z) for z in read_regex_str_qual]
    format_read = partial(_format_read, read_regex_list=read_regex_list,
                          read_template=read_template.read_template,
                          cb_tag=read_template.cb_tag,
                          ub_len=read_template.ub_len,
                          barcode_filter_dict=barcode_filter_dict)

    chunk_size = 8000
    fastq_reader = [read_fastq(fastq_i) for fastq_i in fastq]
    chunks = toolz.partition_all(chunk_size, zip(*fastq_reader))

    num_cpu = mp.cpu_count()
    num_thread = num_thread if num_cpu > num_thread else num_cpu
    seq_chunk_obj = toolz.partition_all(num_thread, chunks)

    fastq_out_all = [fastq_out + str(x) + '.gz' for x in range(num_thread)]
    [gzip.open(x, 'wb').close() for x in fastq_out_all]

    cb_count_all = [cb_count + str(x) + '.csv' for x in range(num_thread)]
    [open(x, 'wt').close() for x in cb_count_all]

    fastq_info = collections.defaultdict(collections.Counter)
    iteration = 0
    results = []

    time_start = time.time()
    pool = mp.Pool(num_thread)
    for fastq_chunk in seq_chunk_obj:
        res = pool.starmap_async(format_read, zip(fastq_chunk, fastq_out_all, cb_count_all))
        results.append(res)

        if len(results) == num_thread * 10:
            results[0].wait()

        while results and results[0].ready():
            iteration += 1
            if not (iteration % 10):
                logger.info(f'Processed {iteration * chunk_size * num_thread:,d} reads!')

            res = results.pop(0)
            chunk_info = res.get()
            _update_fastq_info(fastq_info, chunk_info)

    pool.close()
    pool.join()
    for res in results:
        chunk_info = res.get()
        _update_fastq_info(fastq_info, chunk_info)

    with open('.fastq_count.tsv', 'w') as f:
        for k, v in fastq_info['read'].most_common():
            f.write(f'{k}\t{v}\n')

    cmd_cat_fastq = ' '.join(['cat'] + fastq_out_all + ['>'] + [fastq_out])
    try:
        subprocess.check_output(cmd_cat_fastq, shell=True)
        [os.remove(fastq_file) for fastq_file in fastq_out_all]
    except subprocess.CalledProcessError:
        logger.info(f'Errors in concatenate fastq files')
        sys.exit(-1)
    except OSError:
        logger.info(f'Errors in deleting fastq files')
        sys.exit(-1)

    time_used = time.time() - time_start
    logger.info(f'Formatting fastq done, taking {time_used/3600.0:.3f} hours')
    if not cb_count:
        cb_count = fastq_out + '.cb_count'

    df = _count_cell_barcode_umi(cb_count_all[0])
    for cb_file in cb_count_all[1:]:
        df1 = _count_cell_barcode_umi(cb_file)
        df = pd.concat([df, df1], axis=0)
        df = df.groupby(df.index).sum()
        if df.shape[0] > max_num_cell * 2:
            df = df.sort_values(by=df.columns[0], ascending=False)
            df = df.iloc[:max_num_cell, :]

    try:
        [os.remove(cb_file) for cb_file in cb_count_all]
    except OSError:
        logger.info(f'Errors in deleting cell barcode files')
        sys.exit(-1)

    df = df.sort_values(by=df.columns[0], ascending=False)

    if df.shape[0] > 0:
        df.columns = [str(x) for x in range(df.shape[1])]
        df.index.name = 'cb'
        column_name = list(df.columns.values)
        column_name[0] = 'cb_count'
        df.columns = column_name

        df.to_csv(cb_count, sep='\t')


def _update_fastq_info(fastq_info, chunk_info):
    for fastq_count in chunk_info:
        fastq_info['read'].update(read_pass=fastq_count[0],
                                  read_pass_barcode=fastq_count[1],
                                  read_pass_polyt=fastq_count[2],
                                  read_total=fastq_count[3])


def _count_cell_barcode_umi(cb_file, chunk_size=10 ** 7):
    cb_reader = pd.read_csv(cb_file, header=None, iterator=True,
                            sep='\t', index_col=0)

    chunks = cb_reader.get_chunk(chunk_size)
    chunks = chunks.groupby(chunks.index).sum()

    status = True
    while status:
        try:
            chunk = cb_reader.get_chunk(chunk_size)
            chunks = pd.concat([chunks, chunk], axis=0)
            chunks = chunks.groupby(chunks.index).sum()
        except StopIteration:
            status = False
            logger.info('Read cell barcode counts done.')

    return chunks


def _extract_barcode_pos(barcode_dict, config):
    barcode_reg = []
    pos_all = []

    barcode_filter = dict()
    for barcode_and_pos in barcode_dict:
        barcode, pos = barcode_and_pos
        pos_all.append(pos)

        barcode_reg.append('(?P<' + barcode + '>.{' +
                           str(pos[1] - pos[0] + 1) + '})')
        try:
            value = config[barcode + '_value']
            barcode_filter.update({barcode: ErrorBarcodeHash(value, 1)})
        except KeyError:
            pass

    return barcode_reg, pos_all, barcode_filter


def _extract_input_read_template(read, config):
    read_name = '(@.*)\\n'
    read_plus = '(\\+.*)\\n'
    read_qual = '(.*)\\n'

    filter_dict = dict()
    seq = [(key, value) for key, value in config[read].items()
           if key.startswith('cDNA')]
    if seq:
        read_name = '@(?P<name>.*)\\n'
        read_seq = '(?P<seq>.*)\\n'
        read_qual = '(?P<qual>.*)\\n'
        read_template = read_name + read_seq + read_plus + read_qual
        return read_template, filter_dict, read_template

    cell_barcode = [(key, value) for key, value in config[read].items()
                    if key.startswith('CB') and not key.endswith('value')]

    umi = [(key, value) for key, value in config[read].items()
           if key.startswith('UMI')]

    poly_t = [(key, value) for key, value in config[read].items()
              if key.startswith('polyT')]

    cb_reg, cb_pos, cb_filter = _extract_barcode_pos(cell_barcode, config[read])
    filter_dict.update(cb_filter)

    umi_reg, umi_pos, _ = _extract_barcode_pos(umi, config[read])
    umi_reg = [z.replace('UMI', 'UB') for z in umi_reg]

    pt_reg, pt_pos, _ = _extract_barcode_pos(poly_t, config[read])

    read_pos_start = [z[0] for z in cb_pos]
    read_pos_start += [z[0] for z in umi_pos]
    read_pos_start += [z[0] for z in pt_pos]

    read_pos_end = [z[1] for z in cb_pos]
    read_pos_end += [z[1] for z in umi_pos]
    read_pos_end += [z[1] for z in pt_pos]

    idx = sorted(range(len(read_pos_start)),
                 key=lambda k: read_pos_start[k])

    barcode_tag = cb_reg + umi_reg + pt_reg
    read_pos_start = [read_pos_start[i] for i in idx]
    read_pos_end = [read_pos_end[i] for i in idx]
    barcode_tag = [barcode_tag[i] for i in idx]

    idx_skip = [read_pos_start[i+1] - read_pos_end[i] - 1
                for i in range(0, len(read_pos_start)-1)]

    barcode_skip = ['[ACGTN]{' + str(i) + '}' for i in idx_skip]

    read_seq = barcode_tag[0]
    for i in range(len(read_pos_start)-1):
        if idx_skip[i] == 0:
            read_seq += barcode_tag[i+1]
        else:
            read_seq += barcode_skip[i]
            read_seq += barcode_tag[i+1]

    filter_dict.update(_filter_ploy_t(read_seq))

    if read_pos_start[0] > 1:
        read_seq = '[ACGTN]{' + str(read_pos_start[0]-1) + '}'

    read_seq += '[ACGTN]*'
    read_seq = read_seq + '\\n'

    read_template = read_name + read_seq + read_plus + read_qual

    read_qual = re.sub('>', r'_qual>', read_seq)
    read_qual = re.sub('\[ACGTN\]', '.', read_qual)
    read_template_qual = read_name + read_seq + read_plus + read_qual

    return read_template, filter_dict, read_template_qual


def _filter_ploy_t(read_seq):
    match = re.findall('\?P<polyT>\.{[0-9]+}', read_seq)
    poly_t_count = [int(re.findall(r'\d+', z)[0]) for z in match]

    poly_t_filter = {'polyT': ErrorBarcodeHash('T' * z, 1) for z in poly_t_count}

    return poly_t_filter


def _replace_poly_t(read_seq):
    match = re.findall('\?P<polyT>\.{[0-9]+}', read_seq)

    poly_t_count = [int(re.findall(r'\d+', z)[0]) for z in match]
    poly_t = ['(' + 'T'*z + ')' + '{s<=1}' for z in poly_t_count]

    for z in range(len(match)):
        read_seq = read_seq.replace(match[z], poly_t[z])

    return read_seq


def _infer_read_template(reg_list):
    class ReadInfo(NamedTuple):
        cb: bool
        cb_tag: list
        cb_len: list
        ub: bool
        ub_tag: list
        ub_len: list
        read_template: str

    cb = ub = False
    cb_tag = ub_tag = []
    cb_len = ub_len = []

    read_template = '@'
    reg = ''.join(k for k in reg_list)
    if 'CB' in reg:
        logger.info('Cell barcode in configure file')
        cb = True
        cb_seq_template = _accumulate_barcode('CB', reg)
        cb_template = ':CB_' + cb_seq_template[1]
        read_template += cb_template
        cb_tag = cb_seq_template[0]
        cb_len = cb_seq_template[2]

    if 'UB' in reg:
        logger.info('UMI in config file')
        ub = True
        ub_seq_template = _accumulate_barcode('UB', reg)
        ub_template = ':UB_' + ub_seq_template[1]
        read_template += ub_template
        ub_tag = ub_seq_template[0]
        ub_len = ub_seq_template[2]

    read_template += ':{name}'
    read_template += '\n{seq}\n+\n{qual}\n'

    return ReadInfo(cb=cb, cb_tag=cb_tag, cb_len=cb_len,
                    ub=ub, ub_tag=ub_tag, ub_len=ub_len,
                    read_template=read_template)


def _accumulate_barcode(barcode, seq):
    barcode_num = [sub_str[0] for sub_str in
                   seq.split('?P<' + re.escape(barcode))][1:]
    status = '>' in barcode_num

    barcode_num = ['0' if x == '>' else x for x in barcode_num]
    barcode_num = sorted(barcode_num, key=int)

    if status:
        barcode_num[0] = ''
    barcode_seq = [barcode + num for num in barcode_num]

    barcode_template = ['{' + tag + '}' for tag in barcode_seq]
    barcode_template = '-'.join(barcode_template)

    str_split = 'P<' + barcode + '[0-9]*>.{'
    barcode_len = [sub_str for sub_str in re.split(str_split, seq)][1:]
    barcode_len = [int(re.findall(r'(\d+)', barcode_i)[0])
                   for barcode_i in barcode_len]

    return barcode_seq, barcode_template, barcode_len


def _format_read(chunk, fastq_file, cb_count_file, read_regex_list,
                 read_template, cb_tag, ub_len, barcode_filter_dict):
    reads = []
    num_read = len(chunk)
    num_read_pass = num_read_barcode = num_read_polyt = 0
    num_regex = len(read_regex_list)

    barcode_counter = collections.defaultdict(
        partial(np.zeros, shape=(ub_len[0] + 1), dtype=np.uint32))
    ignore_read = False
    for read_i in chunk:
        read_dict_list = []
        for i, regex_i in enumerate(read_regex_list):
            read_match = regex_i.match(read_i[i])
            if not read_match:
                ignore_read = True
                break

            read_dict_list.append(read_match.groupdict())

        if ignore_read:
            ignore_read = False
            continue

        read1_dict = read_dict_list[0]
        if num_regex > 1:
            for regex_id in range(1, num_regex):
                read1_dict.update(read_dict_list[regex_id])

        cb = [barcode_filter_dict[tag][read1_dict[tag]]
              if tag in barcode_filter_dict.keys() else read1_dict[tag]
              for tag in cb_tag]
        if all(cb):
            cb = '-'.join(cb)
            num_read_barcode += 1
        else:
            ignore_read = True

        ub = read1_dict['UB']
        try:
            poly_t = read1_dict['polyT']
            if not barcode_filter_dict['polyT'][poly_t]:
                ignore_read = True
            else:
                num_read_polyt += 1
        except KeyError:
            pass

        if ignore_read:
            ignore_read = False
            continue
        num_read_pass += 1

        if len(read1_dict['seq']) >= 1:
            read1_dict = read_template.format_map(read1_dict)
            reads.append(read1_dict)

        barcode_counter[cb] += [x == 'T' for x in 'T' + ub]

    with gzip.open(fastq_file, 'ab') as fastq_hd:
        for read in reads:
            fastq_hd.write(bytes(read, 'utf8'))

    df = pd.DataFrame.from_dict(barcode_counter, orient='index')
    if df.shape[0] > 0:
        df = df.sort_values(by=df.columns[0], ascending=False)
        df.index.name = 'cb'
        column_name = list(df.columns.values)
        column_name[0] = 'cb_count'
        df.columns = column_name

        df.to_csv(cb_count_file, sep='\t', mode='a', header=False)

    return num_read_pass, num_read_barcode, num_read_polyt, num_read


def _construct_barcode_regex(bam):
    read_mode = 'r' if bam.endswith('.sam') else 'rb'
    bam_file = AlignmentFile(bam, mode=read_mode)
    first_alignment = next(bam_file)
    bam_file.close()

    barcodes = set()
    for barcode in ['CB_', 'UB_']:
        if barcode in first_alignment.qname:
            barcodes.add(barcode)

    barcode_parser = '.*'
    if 'CB_' in barcodes:
        barcode_parser += ':CB_(?P<CB>[A-Z\-]+)'
    if 'UB_' in barcodes:
        barcode_parser += ':UB_(?P<UB>[A-Z\-]+)'

    if barcode_parser == '.*':
        logger.error('Error: no cell barcodes and UMIs.')
        sys.exit(-1)
    barcode_parser += ':*'

    barcode_parser = re.compile(barcode_parser)
    match = barcode_parser.match(first_alignment.qname)

    cb = _extract_tag(match, 'CB')

    return barcode_parser, cb, read_mode


def _extract_tag(match, tag):
    try:
        tag = match.group(tag)
    except IndexError:
        tag = None

    return tag


def count_feature(*cb, bam, molecular_info_h5, gtf, cb_count, feature_tag='XT:Z',
                  expect_cell=False, force_cell=False, all_cell=False,
                  depth_threshold=1, cell_barcode_whitelist=None):
    """
    Count the number of reads/UMIs mapped to each gene

    :param bam: the input sam/bam file
    :param molecular_info_h5: output the molecular info
    :param cb: the input cell barcode files, can be empty or None
    :param cell_barcode_whitelist: a file contain the selected cell barcodes
    :param gtf: a GTF file
    :param cb_count: a file containing the number of reads mapped to each cell barcode,
           output from format_fastq
    :param feature_tag: the tag representing genes in the input bam file
    :param depth_threshold: only considering UMIs that have at least
           depth_threshold reads support
    :param expect_cell: the expected number of cells in the bam file
    :param force_cell: force to return the number of cells set by expect_cell
    :param all_cell: keep all cell barcodes - can be very slow
    """

    barcode_parser, first_cb, read_mode = _construct_barcode_regex(bam)
    num_cb = len(first_cb.split('-'))
    num_cb_file = len(cb)
    if 0 == num_cb_file:
        cb = [None] * num_cb
    elif num_cb != num_cb_file:
        logger.error(f'Error: the number of input cell barcodes files {num_cb_file} '
                     f'is different from the number of cell barcodes {num_cb} '
                     f'detected in the bam file')

        if num_cb > num_cb_file:
            cb = cb + [None] * (num_cb - num_cb_file)
        else:
            cb = cb[:num_cb]

    # TODO: no cell barcodes detected
    correct_cb_fun, cb_list, cb_remove = _construct_cb_filter(
        cb_count, cb, expect_cell, force_cell, all_cell, cell_barcode_whitelist)

    gene_map_dict = read_gene_map_from_gtf(gtf)

    logger.info('Counting molecular info')
    time_start_count = time.time()

    sam_file = AlignmentFile(bam, mode=read_mode)
    _count_feature_partial = partial(_count_feature,
                                     gene_map_dict=gene_map_dict,
                                     barcode_parser=barcode_parser,
                                     correct_cb_fun=correct_cb_fun,
                                     sam_file=sam_file,
                                     feature_tag=feature_tag)

    track = sam_file.fetch(until_eof=True)
    map_info, read_in_cell, molecular_info = _count_feature_partial(track)

    time_count = time.time() - time_start_count
    logger.info(f'Counting molecular info done - {time_count/3600.0:.3f} hours, '
                f'{int(3600.0 * map_info["num_alignment"]/time_count):,d} '
                f'alignments/hour\n')

    # TODO: still output results
    if len(molecular_info) == 0:
        logger.error('Error: no reads mapped to features.')
        sys.exit(-1)

    name = ['cell',
            'gene',
            'umi',
            'depth',
            ]

    logger.info('Converting to a dataframe')
    convert_time = time.time()
    molecular_info = pd.Series(molecular_info).reset_index()
    molecular_info.columns = name

    for col in name[:3]:
        molecular_info.loc[:, col] = molecular_info[col].astype('category')

    convert_time = time.time() - convert_time
    logger.info(f'Converting to a dataframe done, '
                f'taking {convert_time/60.0:.3f} minutes\n')

    molecular_info.columns = name
    if num_cb > 1 and cb_list:
        molecular_info = molecular_info.loc[molecular_info['cell'].isin(cb_list), :]

    if cb_remove:
        molecular_info = molecular_info.loc[~molecular_info['cell'].isin(cb_remove), :]

    molecular_info = molecular_info.loc[molecular_info['depth'] >= 0.95, :]
    molecular_info['depth'] = \
        np.floor(molecular_info['depth'].values + 0.5).astype('uint32')

    molecular_info = molecular_info.sort_values(name[:3])
    molecular_info = molecular_info.reset_index(drop=True)

    map_info = pd.Series(map_info)
    read_in_cell = pd.DataFrame.from_dict(read_in_cell, orient='index')

    logger.info('Writing molecular info')
    write_time = time.time()

    feature = gene_map_dict.values()
    feature = pd.Series(index=set(feature))
    feature = feature.sort_index()

    with pd.HDFStore(molecular_info_h5, mode='w') as hf:
        hf.put('molecular_info', molecular_info, format='table', data_columns=True)
        hf.put('map_info', map_info)
        hf.put('feature', feature)
        hf.put('read_in_cell', read_in_cell)
    del molecular_info

    write_time = time.time() - write_time
    logger.info(f'Writings molecular info done, '
                f'taking {write_time/60.0:.3f} minutes\n')

    _convert_count_to_matrix(molecular_info_h5, molecular_info_h5,
                             depth_threshold=depth_threshold)


def _count_feature(track, gene_map_dict, barcode_parser,
                   correct_cb_fun, sam_file, feature_tag='XT:Z'):
    search_undetermined = re.compile('N').search

    read_name = None
    feature_tag_value_pre = None
    filt_multiple_gene_barcode = False
    count_read = False
    cb_i = feature_tag_value = ub_i = None
    num_aln_read = 0
    pass_filter = False

    map_info = collections.defaultdict(int)
    read_in_cell = collections.Counter()
    molecular_info = collections.defaultdict(int)

    for aln in track:
        if map_info['num_alignment'] and not map_info['num_alignment'] % 10000000:
            logger.info(f'Parsed {map_info["num_alignment"]:,d} alignments.')
            logger.info(f'{map_info["num_unique_read"]:,d} unique reads, '
                        f'{map_info["num_count_read"]:,d} reads kept.')

            logger.info(f'{map_info["num_unmapped_read"]:,d} unmapped reads were filtered.')
            logger.info(f'{map_info["num_barcode_with_na"]:,d} reads '
                        f'were filtered for including NA in barcodes.\n')

        num_aln_read_pre = num_aln_read

        filter_read_unmapped = False
        filter_read_na = False
        filter_read_barcode = False
        map_info['num_alignment'] += 1

        num_aln_read = aln.get_tag('NH')
        new_read = aln.qname != read_name
        if new_read:
            read_name = aln.qname
            if count_read:
                map_info['num_count_read'] += 1
                record_tuple = (cb_i, feature_tag_value, ub_i)
                molecular_info[record_tuple] += 1
            elif pass_filter and (num_aln_read_pre > 1):
                map_info['num_barcode_with_na'] += 1

            pass_filter = False
            count_read = False
            filt_multiple_gene_barcode = True
            feature_tag_value_pre = None

            map_info['num_unique_read'] += 1
            if num_aln_read == 0:
                map_info['num_unmapped_read'] += 1
                filter_read_unmapped = True

            # check cb
            match = barcode_parser.match(aln.qname)
            cb_i = _extract_tag(match, 'CB')
            cb_i_list = cb_i.split('-')
            num_na_in_cb = _count_not_specified(cb_i_list)
            if any(num_na_in_cb > 1) or sum(num_na_in_cb) > len(num_na_in_cb):
                filter_read_na = True

            cb_i = correct_cb_fun(cb_i_list)
            if cb_i:
                read_in_cell[cb_i] += 1
            elif not aln.is_unmapped:
                map_info['num_barcode_with_na'] += 1
                filter_read_barcode = True

            if filter_read_unmapped or filter_read_na or filter_read_barcode:
                count_read = False
                continue

            ub_i = _extract_tag(match, 'UB')
            if ub_i and search_undetermined(ub_i):
                map_info['num_barcode_with_na'] += 1
                continue

            if ub_i and ub_i == len(ub_i) * ub_i[0]:
                map_info['num_barcode_with_na'] += 1
                continue

            try:
                feature_tag_value = aln.get_tag(feature_tag)
            except KeyError:
                feature_tag_value = sam_file.getrname(aln.reference_id)
                if aln.get_tag('XS:Z') == 'Unassigned_Ambiguity':
                    map_info['num_barcode_with_na'] += 1
                    continue

            pass_filter = True
            filt_multiple_gene_barcode = False

            try:
                feature_tag_value = gene_map_dict[feature_tag_value]
            except KeyError:
                if num_aln_read == 1:
                    map_info['num_barcode_with_na'] += 1
                continue

            feature_tag_value_pre = feature_tag_value

            count_read = True
        else:
            if filt_multiple_gene_barcode:
                continue

            try:
                feature_tag_value = aln.get_tag(feature_tag)
            except KeyError:
                feature_tag_value = sam_file.getrname(aln.reference_id)
                if aln.get_tag('XS:Z') == 'Unassigned_Ambiguity':
                    filt_multiple_gene_barcode = True
                    count_read = False
                    continue

            try:
                feature_tag_value = gene_map_dict[feature_tag_value]
            except KeyError:
                feature_tag_value = feature_tag_value_pre
                continue

            if feature_tag_value_pre and feature_tag_value_pre != feature_tag_value:
                filt_multiple_gene_barcode = True
                count_read = False
                continue
            feature_tag_value_pre = feature_tag_value

            count_read = True  # with valid feature_tag_value

    if count_read:
        map_info['num_count_read'] += 1
        record_tuple = (cb_i, feature_tag_value, ub_i)
        molecular_info[record_tuple] += 1

    return map_info, read_in_cell, molecular_info


def _construct_cb_filter(cb_count, cb, expect_cell, force_cell,
                         all_cell, cell_barcode_whitelist):
    cb_list = []
    cb_remove = []
    if all_cell:
        correct_cb_fun = _filter_tag_fun(cb, max_distance=1, correct=True)
    else:
        if cell_barcode_whitelist:
            with open(cell_barcode_whitelist, 'r') as file_handle:
                cb_list = [line.rstrip('\n') for line in file_handle]
        else:
            cb_list, cb_remove = _get_candidate_barcode(cb_count, cb,
                                                        expect_cell=expect_cell,
                                                        force_cell=force_cell)

        num_cell = len(cb_list)
        logger.info(f'Detected {num_cell:,d} candidate cell barcodes.')
        if num_cell <= 0:
            sys.exit(-1)

        cb_list_split = [cb.split('-') for cb in cb_list]
        cb_df = pd.DataFrame(cb_list_split)

        cb_list_split = [''] * len(cb_df.columns)
        for cb in cb_df:
            cb_list_split[cb] = cb_df[cb].unique()

        if len(cb_df.columns) > 1:
            barcode_hash = _create_barcode_hash(cb_list)
            cb_hash = [ErrorBarcodeHashConstraint(cb, barcode_hash[idx])
                       for idx, cb in enumerate(cb_list_split)]
            correct_cb_fun = partial(_filter_tag_multi, tag_hash=cb_hash, correct=True)
        else:
            cb_hash = [ErrorBarcodeHash(cb) for cb in cb_list_split]
            correct_cb_fun = partial(_filter_tag, tag_hash=cb_hash, correct=True)

    return correct_cb_fun, cb_list, cb_remove


def _count_not_specified(barcode):
    if not barcode:
        return np.array([0])

    na_count = [barcode_i.count('N') for barcode_i in barcode]

    return np.array(na_count)


def _get_candidate_barcode(cb_count_file, cb_file,
                           plot_prefix='.',
                           expect_cell=False,
                           force_cell=False):
    cb = pd.read_csv(cb_count_file, sep='\t')
    cb_name = cb['cb'].str.split('-', expand=True)

    cb_len = [len(z) for z in cb_name.iloc[0, :]]
    num_cb = len(cb_len)

    idx = False
    for cb_idx in range(num_cb):
        filt_cb = [cb_char * cb_len[cb_idx] for cb_char in ['N', 'G']]
        idx = cb_name.iloc[:, 0].isin(filt_cb) | idx
    cb = cb.loc[~idx, :]
    cb = cb.reset_index(drop=True)

    cb_count = dict(zip(cb.cb, cb.cb_count))

    candidate_cb_whitelist = get_cell_whitelist(cb_count,
                                                plot_prefix=plot_prefix,
                                                expect_cell=expect_cell,
                                                force_cell=force_cell)
    candidate_cb_whitelist_refine = \
        _refine_whitelist(list(candidate_cb_whitelist), cb_file)

    merge_barcode = [x != 'None' and x for x in cb_file]
    if any(merge_barcode):
        merge_barcode = False

    cb_list, cb_remove = _merge_del_barcode(candidate_cb_whitelist_refine,
                                            barcode_count=cb, min_distance=1,
                                            merge_barcode=merge_barcode)

    with open('._detected_cb.tsv', 'wt') as file_handle:
        for cb_whitelist in np.setdiff1d(cb_list, cb_remove):
            file_handle.write(f'{cb_whitelist}\n')

    return cb_list, cb_remove


def _refine_whitelist(cb_whitelist, cb_file=None, max_na_per_cb=1):
    cb_hash = []
    if cb_file is not None:
        cb_file = list(cb_file)
        cb_file = [None if cb_i == 'None' else cb_i for cb_i in cb_file]
        cb_hash = _construct_hash(cb_whitelist, cb_file)

    num_cb = len(cb_whitelist[0].split('-'))
    cb_whitelist_corrected = collections.defaultdict(set)
    for cell_barcode in list(cb_whitelist):
        cell_barcode_list = cell_barcode.split('-')
        na_count = _count_not_specified(cell_barcode_list)
        if any(na_count > max_na_per_cb) or sum(na_count) > num_cb:
            continue

        cb = cell_barcode
        if any(cb_hash):
            cb = _correct_cell_barcode(cell_barcode_list, cb_hash)
            if any(cb):
                cb = '-'.join(cb)
            else:
                continue

        cb_whitelist_corrected[cell_barcode].add(cb)

    return cb_whitelist_corrected


def _construct_hash(cb_whitelist, tag_file):
    num_tag = len(tag_file)
    tag_hash = [''] * num_tag

    # Add comments
    for i in range(num_tag):
        if tag_file[i]:
            with open(tag_file[i], 'r') as file_handle:
                tag_i = [line.rstrip('\n') for line in file_handle]

            if len(tag_i) > 5000:
                cell_barcode_list = []
                for cell_barcode in list(cb_whitelist):
                    cell_barcode_list.append(cell_barcode.split('-')[i])

                white_list_map = \
                    _generate_barcode_whitelist_map(cell_barcode_list, tag_i, 1)
                tag_i = [list(v)[0] for k, v in white_list_map.items()]
                tag_i = list(set(tag_i))

            tag_hash[i] = ErrorBarcodeHash(tag_i, edit_distance=1)

    return tag_hash


def _correct_cell_barcode(cell_barcode, cb_hash):
    num_cb = len(cb_hash)

    cb_corrected = cell_barcode
    for i in range(num_cb):
        if cb_hash[i]:
            candidate_cb = cb_hash[i][cell_barcode[i]]
            if candidate_cb:
                cb_corrected[i] = candidate_cb
            else:
                return [None]

    return cb_corrected


def _generate_barcode_whitelist_map(barcode, whitelist, min_distance=1):
    barcode_to_whitelist = collections.defaultdict(set)

    whitelist = set([str(x).encode('utf-8') for x in whitelist])

    num_cpu = mp.cpu_count()
    pool = mp.Pool(num_cpu)

    _partial_map_single_barcode_to_whitelist = \
        partial(_map_single_barcode_to_whitelist, whitelist=whitelist,
                min_distance=min_distance)

    corrected_barcode = pool.map(_partial_map_single_barcode_to_whitelist, barcode)

    for idx, barcode_i in enumerate(corrected_barcode):
        if barcode_i is not None:
            barcode_to_whitelist[barcode[idx]].add(barcode_i)

    return barcode_to_whitelist


def _map_single_barcode_to_whitelist(barcode, whitelist, min_distance=1):
    match = None
    barcode_in_bytes = str(barcode).encode('utf-8')

    for white_barcode in whitelist:
        if barcode_in_bytes in whitelist:
            match = barcode
            break

        if compute_edit_distance(barcode_in_bytes, white_barcode) <= min_distance:
            if match is not None:
                logging.info(f'Warning: barcode {str(barcode)} can be '
                             f'mapped to more than one candidate barcodes')
                match = None
                break
            else:
                match = white_barcode.decode('utf-8')

    return match


def _merge_del_barcode(barcode_dict, barcode_count, min_distance=1, merge_barcode=False):
    barcode_list = list(barcode_dict.keys())

    idx = barcode_count.cb.isin(barcode_list)
    barcode_count_filt = barcode_count.loc[idx, :]

    barcode_corr = [barcode_dict[x] for x in barcode_count_filt.cb]
    idx = [len(x) > 0 for x in barcode_corr]

    barcode_count_filt = barcode_count_filt.iloc[idx, :]

    barcode_corr = list(chain(*barcode_corr))
    barcode_count_filt.cb = barcode_corr

    barcode_count_filt = barcode_count_filt.groupby('cb').sum()

    umi_len = barcode_count_filt.shape[1]
    barcode_count_filt_ratio = barcode_count_filt.iloc[:, 1:umi_len].div(
        barcode_count_filt.cb_count, axis=0)

    idx = barcode_count_filt_ratio.gt(0.80, axis=0)
    idx = idx | barcode_count_filt_ratio.lt(0.005, axis=0)
    count_indel = idx.sum(axis=1)

    if sum(count_indel == 0) <= 100000 and merge_barcode:
        barcode_whitelist = \
            _merge_corrected_barcode(barcode_count_filt.loc[count_indel == 0, :])
    else:
        barcode_whitelist = barcode_count_filt_ratio.index[count_indel == 0].tolist()

    barcode_correctable = barcode_count_filt_ratio.index[count_indel == 1].tolist()

    whitelist_remove = []
    if len(barcode_correctable) > 0:
        barcode_corrected = _correct_del_barcode(
            barcode_count_filt.loc[barcode_correctable, :], min_distance)

        barcode_corrected_list = list(barcode_corrected.keys())
        barcode_corrected_list_mut = [x[:-1]+'N' for x in barcode_corrected_list]
        whitelist_dist = [_find_neighbour_barcode(x, barcode_corrected_list_mut, 1)
                          for x in barcode_whitelist]

        whitelist_remove = [barcode_whitelist[k] for k, v in enumerate(whitelist_dist)
                            if len(v[0]) > 0]

        barcode_whitelist.extend(barcode_corrected_list)

    return barcode_whitelist, whitelist_remove


# N**2 complexity
def _merge_corrected_barcode(barcode_count):
    barcode_count = barcode_count.sort_values('cb_count', ascending=False)

    barcode = barcode_count.index.astype(str)
    barcode_coverage = dict(zip(barcode, barcode_count.cb_count))

    barcode_list = collections.deque()
    barcode_list.append(barcode[0])

    num_barcode = len(barcode)
    if num_barcode <= 1:
        return barcode_list

    for barcode_i in barcode[1:]:
        idx = _find_neighbour_barcode(barcode_i, barcode_list)

        num_neighbour = len(idx[0])
        if num_neighbour == 0:
            barcode_list.append(barcode_i)
            continue
        elif num_neighbour == 1:
            candidate_barcode_idx = idx[0][0]
            candidate_barcode_coverage = \
                barcode_coverage[barcode_list[candidate_barcode_idx]]

            if barcode_coverage[barcode_i] > candidate_barcode_coverage / 10.0:
                barcode_list.append(barcode_i)
                continue

    return barcode_list


def _find_neighbour_barcode(barcode, barcode_list, min_distance=1):
    edit_dist = np.array([compute_edit_distance(barcode.encode('utf-8'), x.encode('utf-8'))
                         for x in barcode_list])
    idx_filt = np.where(edit_dist <= min_distance)

    return idx_filt


def _correct_del_barcode(barcode_count, min_distance=1):
    barcode_count = barcode_count.sort_values('cb_count', ascending=False)
    barcode_all = np.asarray(barcode_count.index.tolist())

    barcode_whitelist = collections.defaultdict(set)
    while len(barcode_all):
        barcode_i = barcode_all[0]
        barcode_whitelist[barcode_i].add(barcode_i)
        barcode_all = barcode_all[1:]

        idx = _find_neighbour_barcode(barcode_i, barcode_all, min_distance)
        if len(idx[0]):
            barcode_ = barcode_all[idx]
            barcode_whitelist[barcode_i].update(list(barcode_))
            barcode_all = np.delete(barcode_all, idx)

    return barcode_whitelist


def _create_barcode_hash(barcode):
    barcode_split = [barcode_i.split('-') for barcode_i in barcode]

    barcode_df = pd.DataFrame(barcode_split, barcode)
    num_barcode = len(barcode_split[0])

    barcode_split_uniq = [''] * len(barcode_df.columns)
    for barcode_i in barcode_df:
        barcode_split_uniq[barcode_i] = barcode_df[barcode_i].unique()

    barcode_hash = [collections.defaultdict(list) for _ in range(num_barcode)]
    for i in range(num_barcode):
        barcode_i = barcode_split_uniq[i]
        for barcode_ii in barcode_i:
            idx = barcode_df[i] == barcode_ii
            barcode_hash[i][barcode_ii] = set(barcode_df.index[idx].tolist())

    return barcode_hash


def convert_count_to_matrix(molecular_info, out_prefix, depth_threshold):
    _convert_count_to_matrix(molecular_info, out_prefix, depth_threshold)


def _convert_count_to_matrix(molecular_info, out_prefix, depth_threshold):
    feature = pd.read_hdf(molecular_info, key='feature')
    molecular_info = pd.read_hdf(molecular_info, key='molecular_info')

    logger.info('Collapsing UMIs')
    write_time = time.time()
    molecular_info = _collapse_umi(molecular_info)
    write_time = time.time() - write_time
    logger.info(f'Collapsing UMIs done, taking {write_time/60.0:.3f} minutes')

    df = _generate_fake_count(molecular_info.iloc[0, :],
                              feature, depth=depth_threshold+1)
    molecular_info = pd.concat([df, molecular_info], ignore_index=True)

    num_gene = len(feature)
    _transform_write_sparse_matrix(molecular_info, num_gene,
                                   sum_type='umi', out_prefix=out_prefix,
                                   depth_threshold=depth_threshold)
    _transform_write_sparse_matrix(molecular_info, num_gene,
                                   sum_type='transcript', out_prefix=out_prefix,
                                   depth_threshold=depth_threshold)


def _transform_write_sparse_matrix(molecular_info, num_gene,
                                   sum_type, out_prefix, depth_threshold):

    logger.info('Converting to sparse matrix')
    convert_time = time.time()

    query_filer = f'depth >= {depth_threshold}'
    if 'umi' == sum_type:
        base_name = out_prefix + '_read'
        count_collapsed = molecular_info.groupby(
            ['cell', 'gene'])
        count_collapsed = count_collapsed['depth'].sum()

        count_collapsed[:num_gene] -= (depth_threshold + 1)
        count_collapsed += 0.5
        count_collapsed = count_collapsed.astype(int)
    else:
        base_name = out_prefix + '_depth_' + str(depth_threshold) + '_transcript'
        count_collapsed = molecular_info.query(query_filer).groupby(
            ['cell', 'gene'])
        count_collapsed = count_collapsed['umi'].size()

        count_collapsed[:num_gene] -= 1

    del molecular_info

    count, count_row_name, count_column_name = _convert_to_coo(count_collapsed)
    del count_collapsed

    convert_time = time.time() - convert_time
    logger.info(f'Converting to a sparse matrix done, '
                f'taking {convert_time/60.0:.3f} minutes')

    logger.info('Output results')
    write_time = time.time()
    pd.Series(count_row_name).to_csv(base_name + '_gene.tsv',
                                     index=False, header=False)
    pd.Series(count_column_name).to_csv(base_name + '_barcode.tsv',
                                        index=False, header=False)

    if 'umi' == sum_type:
        with open(base_name + '.mtx', 'w+b') as out_handle:
            scipy.io.mmwrite(out_handle, count)
    else:
        with open(base_name + '.mtx', 'w+b') as out_handle:
            scipy.io.mmwrite(out_handle, count)

    write_time = time.time() - write_time
    logger.info(f'Writing final results done, '
                f'taking {write_time/60.0:.3f} minutes')


def _map_barcode_to_whitelist(barcode, whitelist, min_distance=1):
    whitelist = set([str(x).encode('utf-8') for x in whitelist])

    iter_i = 0
    for barcode_i in barcode:
        match = barcode_i
        barcode_in_bytes = str(barcode_i).encode('utf-8')

        for white_barcode in whitelist:
            if barcode_in_bytes in whitelist:
                break

            if compute_edit_distance(barcode_in_bytes, white_barcode) <= min_distance:
                match = white_barcode.decode('utf-8')
                break

        barcode[iter_i] = match
        iter_i += 1

    return barcode


def _collapse_barcode_edit(barcode, value, min_distance=1):
    id_srt = value.argsort()[::-1]

    barcode = barcode[id_srt]
    value = value[id_srt]

    max_barcode = value[0]
    threshold = max_barcode * 0.1
    if threshold < 2:
        threshold = 2
    elif threshold > 5:
        threshold = 5
    id_whitelist = value > threshold

    whitelist_candidate = barcode[id_whitelist]
    noise_candidate = barcode[~id_whitelist]

    if len(noise_candidate) > 0 and len(whitelist_candidate) > 0:
        corrected_noise = _map_barcode_to_whitelist(noise_candidate,
                                                    whitelist_candidate,
                                                    min_distance)

        barcode[~id_whitelist] = corrected_noise

    return barcode, value


def _collapse_umi(x, min_distance=1):
    id_start = x.duplicated(['cell', 'gene'])
    id_start = id_start[id_start == False].index.tolist()

    id_end = id_start[1:]
    id_end.append(x.shape[0])

    value = x['depth'].values
    umi = x['umi'].values.astype('str')

    for gene in np.arange(len(id_end)):
        id_gene = np.arange(id_start[gene], id_end[gene])

        if len(id_gene) <= 1:
            continue

        umi_gene = umi[id_gene]
        value_gene = value[id_gene]

        umi_gene, _ = _collapse_barcode_edit(umi_gene, value_gene, min_distance)
        umi[id_gene] = umi_gene

    x['umi'] = umi
    x = x.groupby(['cell', 'gene', 'umi'])['depth'].sum()
    x = x.reset_index(drop=False)

    return x


def _convert_to_coo(data_series):
    data_sp = data_series.to_sparse()

    data_sp, row_name, column_name = data_sp.to_coo(column_levels=['cell'],
                                                    row_levels=['gene'])
    coo_tuple = collections.namedtuple('coo_tuple', ['x', 'row_name', 'column_name'])

    return coo_tuple(data_sp, row_name, column_name)


def _generate_fake_count(row_of_df, feature, depth=1.5):
    index_name = row_of_df.index
    df = pd.DataFrame(columns=index_name)
    df['gene'] = feature.index.values
    df['umi'] = 'N' * len(row_of_df['umi'])
    df['depth'] = depth

    for index_name_other in list(set(index_name) - {'gene', 'umi', 'depth'}):
        df[index_name_other] = row_of_df[index_name_other]

    return df


def down_sample(molecular_info, total_read=None, total_cell=None, mean_read=None,
                out_prefix='.', depth_threshold=1, seed=0):
    """
    Down-sampling the molecular_info such that each library has the same number of reads

    :param molecular_info: molecular_info: the input molecular info data frame
    :param total_read: the total number of reads for these libraries
    :param total_cell: the total number of cells
    :param mean_read: the expected number of reads per cell after down-sampling
    :param out_prefix: the prefix of the output matrices
    :param depth_threshold: the coverage threshold to consider
    :param seed used for random sampling
    """

    feature = pd.read_hdf(molecular_info, key='feature')
    output_molecular_info = pd.read_hdf(molecular_info, key='molecular_info')
    for col in output_molecular_info.columns[:3]:
        output_molecular_info.loc[:, col] = output_molecular_info[col].astype('category')

    output_molecular_info = output_molecular_info.loc[
                            output_molecular_info['depth'] >= 1, :]
    output_molecular_info.reset_index(drop=True, inplace=True)

    map_info = pd.read_hdf(molecular_info, key='map_info')
    if total_read is None:
        total_read = map_info['num_unique_read']
    else:
        total_read = max(total_read, map_info['num_unique_read'])

    read_in_cell = pd.read_hdf(molecular_info, key='read_in_cell')
    if total_cell is None:
        total_cell = read_in_cell.shape[0]
    else:
        total_cell = min(total_cell, read_in_cell.shape[0])

    if mean_read is None:
        mean_read = (10000, )
    elif not isinstance(mean_read, tuple):
        mean_read = (mean_read, )

    cell_vec = output_molecular_info['depth'].copy()
    for mean_read_i in mean_read:
        random.seed(seed)
        seed += 1

        _down_sample(output_molecular_info, feature, total_read, total_cell,
                     mean_read_i, out_prefix, depth_threshold=depth_threshold)
        output_molecular_info['depth'] = cell_vec


def _down_sample(molecular_info, feature, total_read, total_cell, mean_read,
                 out_prefix, depth_threshold=1):
    expect_read = mean_read * total_cell
    if expect_read > total_read:
        return ()

    cell_vec = molecular_info['depth']
    value = cell_vec.tolist()
    id_end = np.cumsum(value)
    id_start = np.append(0, id_end)
    id_start = id_start[:-1]

    read_num_subsample = (id_end[-1] / total_read * 1.0) * expect_read
    read_num_subsample = int(read_num_subsample + 0.5)

    id_keep = sorted(random.sample(range(id_end[-1]), read_num_subsample))
    expanded_count = np.zeros(id_end[-1], dtype=np.int32)
    expanded_count[id_keep] = 1

    value = _add_umis(expanded_count, id_start, id_end)
    molecular_info['depth'] = value
    output_molecular_info = molecular_info.loc[molecular_info['depth'] >= 1, :].copy()
    output_molecular_info.reset_index(drop=True, inplace=True)

    logger.info('Collapsing UMIs')
    write_time = time.time()
    output_molecular_info = _collapse_umi(output_molecular_info)
    write_time = time.time() - write_time
    logger.info(f'Collapsing UMIs done, taking {write_time / 60.0:.3f} minutes')

    out_prefix = out_prefix + '_sample_' + str(mean_read)
    _calculate_cell_gene_matrix(output_molecular_info, feature,
                                out_prefix=out_prefix,
                                depth_threshold=depth_threshold)


def down_sample_cell(molecular_info, expect_read=None, out_prefix='', depth_threshold=1):
    """
    Down-sampling the molecular_info such that each cell has the same number of reads

    :param molecular_info: the input molecular info data frame
    :param expect_read: each cell to have expect_read
    :param out_prefix: the prefix of the output matrices
    :param depth_threshold: the coverage threshold to consider
    """

    feature = pd.read_hdf(molecular_info, key='feature')
    output_molecular_info = pd.read_hdf(molecular_info, key='molecular_info')

    for col in output_molecular_info.columns[:3]:
        output_molecular_info.loc[:, col] = output_molecular_info[col].astype('category')

    output_molecular_info = output_molecular_info.loc[
                            output_molecular_info['depth'] >= 1, :]
    name = output_molecular_info.columns.tolist()
    name = name[:-1]
    output_molecular_info = output_molecular_info.sort_values(name)  # already sorted?

    read_in_cell = pd.read_hdf(molecular_info, key='read_in_cell')

    if expect_read is None:
        expect_read = (10000, )
    elif not isinstance(expect_read, tuple):
        expect_read = (expect_read, )
    for num_read in expect_read:
        _down_sample_cell(output_molecular_info, feature, read_in_cell, num_read,
                          out_prefix, depth_threshold=depth_threshold)


def _down_sample_cell(molecular_info, feature, read_in_cell, expect_read,
                      out_prefix, depth_threshold=1):

    read_in_cell = read_in_cell[read_in_cell[0] >= expect_read]

    num_cell = read_in_cell.shape[0]
    if num_cell == 0:
        logger.error(f'All cells have less than {expect_read:,d} reads, aborting.')
        sys.exit(-1)

    cell = set(read_in_cell.index.tolist())
    output_molecular_info = molecular_info.loc[
        molecular_info['cell'].isin(cell)].copy()

    # Update the values of output_molecular_info
    output_molecular_info.reset_index(drop=True, inplace=True)

    cell_start = output_molecular_info.drop_duplicates('cell')
    cell = cell_start['cell'].tolist()
    cell_start = cell_start.index.tolist()

    cell_end = cell_start[1:]
    cell_end.append(output_molecular_info.shape[0])

    read_in_cell = read_in_cell.to_dict('dict')[0]

    umi_count = output_molecular_info['depth'].as_matrix()
    umi_count = umi_count.reshape([umi_count.shape[0], 1])
    random.seed(0)

    for idx_i, cell_i in enumerate(cell):
        logger.info(f'Subsample cell {idx_i:,d}')

        id_case, id_end = _find_pos(cell_start[idx_i], cell_end[idx_i], umi_count)
        id_start = np.append(0, id_end)
        id_start = id_start[:-1]

        num_umi = int(expect_read * id_end[-1] / read_in_cell[cell_i] + 0.5)
        if id_end[-1] < num_umi:
            continue

        id_keep = sorted(random.sample(range(id_end[-1]), num_umi))
        expanded_count = np.zeros(id_end[-1], dtype=np.int32)
        expanded_count[id_keep] = 1

        value = _add_umis(expanded_count, id_start, id_end)
        umi_count[id_case] = value

    logger.info('Collapsing UMIs')
    write_time = time.time()
    output_molecular_info = _collapse_umi(output_molecular_info)
    write_time = time.time() - write_time
    logger.info(f'Collapsing UMIs done, taking {write_time/60.0:.3f} minutes')

    out_prefix = out_prefix + '.sample_' + str(expect_read)
    _calculate_cell_gene_matrix(output_molecular_info, feature, out_prefix=out_prefix,
                                depth_threshold=depth_threshold)


@numba.jit(nopython=True)
def _add_umis(x, id_start, id_end):
    y = np.zeros((len(id_start), 1), dtype=np.uint32)
    for i in range(len(id_start)):
        y[i] = np.sum(x[id_start[i]:id_end[i]])

    return y


@numba.jit(nopython=True)
def _find_pos(cell_start_index, cell_end_index, umi_count):
    id_case = np.arange(cell_start_index, cell_end_index)

    value = umi_count[id_case]
    id_end = np.cumsum(value)

    return id_case, id_end


def _calculate_cell_gene_matrix(molecular_info, feature, out_prefix, depth_threshold):
    df = _generate_fake_count(molecular_info.iloc[0, :],
                              feature, depth=depth_threshold+1)
    molecular_info = pd.concat([df, molecular_info], ignore_index=True)

    num_gene = len(feature)
    _transform_write_sparse_matrix(molecular_info, num_gene,
                                   sum_type='transcript', out_prefix=out_prefix,
                                   depth_threshold=depth_threshold)

    _transform_write_sparse_matrix(molecular_info, num_gene,
                                   sum_type='umi', out_prefix=out_prefix,
                                   depth_threshold=depth_threshold)


def _filter_tag_fun(tag_file, max_distance, correct=True):
    if tag_file is not None:
        tag_file = list(tag_file)
        tag_file = [None if x == 'None' else x for x in tag_file]
    else:
        tag_file = [None]

    tag_hash = [None] * len(tag_file)
    for i, tag_file_i in enumerate(tag_file):
        if not tag_file_i:
            continue

        with open(tag_file_i, 'r') as file_handle:
            tag_seq = {line.rstrip('\n') for line in file_handle}
        tag_hash[i] = ErrorBarcodeHash(tag_seq, max_distance)

    filter_tag_fun = partial(_filter_tag, tag_hash=tag_hash, correct=correct)

    return filter_tag_fun


def _filter_tag(tag, tag_hash, correct=True):
    tag_corrected = list(tag)

    for i, tag_i in enumerate(tag):
        if not tag_hash[i]:
            continue

        tag_corrected[i] = tag_hash[i][tag_i]
        if not tag_corrected[i]:
            return None

    if correct:
        return '-'.join(tag_corrected)
    else:
        return '-'.join(tag)


def _filter_tag_multi(tag, tag_hash, correct=True):
    tag_corrected = list(tag)

    tag_full = '-'.join(tag)
    for i, tag_i in enumerate(tag):
        if not tag_hash[i]:
            continue

        tag_corrected[i] = tag_hash[i][tag_i, tag_full]
        if not tag_corrected[i]:
            return None

    if correct:
        return '-'.join(tag_corrected)
    else:
        return '-'.join(tag)


def annotate_bam(bam, gtf, featureCounts='featureCounts',
                 annotate_multi_mapping=True, annotate_intron=False,
                 strand=1, num_thread=4):
    import subprocess

    logger.info('Running featureCounts to annotate alignments.')
    start_time = time.time()

    cmd_feature_count = ' '.join([featureCounts,
                                  '-g gene_id', '-t exon', '-R BAM',
                                  '-T', str(num_thread),
                                  '-F GTF', '-a', gtf, '-s', str(strand),
                                  ])

    if annotate_multi_mapping:
        cmd_feature_count = ' '.join([cmd_feature_count, '-M'])

    try:
        cmd = ' '.join([cmd_feature_count, '-o', bam + '.annotation', bam])
        print(cmd, '\n')

        subprocess.check_output(cmd, shell=True)
    except subprocess.CalledProcessError:
        logger.info(f'Running featureCount errors')
        sys.exit(-1)
    except OSError:
        logger.info(f'featureCounts not found or cannot run!'
                    f' Please specify which featureCounts to use.')
        sys.exit(-1)

    if annotate_intron:
        cmd_feature_count = re.sub('-t exon', '-t transcript', cmd_feature_count)
        cmd = ' '.join([cmd_feature_count, '-o',
                        bam + '.annotation.intron',
                        bam + '.featureCounts.bam'])
        try:
            subprocess.check_output(cmd, shell=True)

            cmd = ' '.join(['mv', bam + '.featureCounts.bam' + '.featureCounts.bam',
                            bam + '.featureCounts.bam'])
            subprocess.check_output(cmd, shell=True)
        except subprocess.CalledProcessError:
            logger.info(f'Running featureCount errors')
            sys.exit(-1)
        except OSError:
            logger.info(f'featureCounts not found or cannot run! ' 
                        f'Please specify which featureCounts to use.')
            sys.exit(-1)

    feature_count_time = time.time() - start_time
    logger.info(f'Annotating features done successfully, '
                f'taking {feature_count_time/60.0:.3f} minutes')
