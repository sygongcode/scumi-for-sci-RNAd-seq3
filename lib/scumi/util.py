import numba
import logging
import regex as re
import gzip

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s: %(levelname)s: %(message)s')
logger = logging.getLogger(__name__)


@numba.jit(nopython=True)
def compute_edit_distance(str1, str2):
    edit_distance = 0
    for char1, char2 in zip(str1, str2):
        if char1 != char2:
            edit_distance += 1

    return edit_distance


def read_gene_map_from_gtf(gtf):
    logger.info('Read GTF {}'.format(gtf))

    gtf_open = gzip.open if gtf.endswith('.gz') else open

    gene_map = {}
    with gtf_open(gtf, 'r') as gtf_file:
        for line in gtf_file:
            if line.startswith('#'):
                continue
            field = line.strip().split('\t')

            if field[2] != 'gene':
                continue

            try:
                gene_id = re.search('gene_id "(.+?)"', field[-1]).group(1)
                gene_name = re.search('gene_name "(.+?)"', field[-1]).group(1)
            except AttributeError:
                continue

            gene_map[gene_id] = gene_id + '_' + gene_name

    logger.info(f'Parse GTF {gtf} done')

    return gene_map
