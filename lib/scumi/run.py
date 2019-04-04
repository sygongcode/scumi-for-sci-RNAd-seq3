from . import scumi


def merge_fastq(args):
    config = args.config
    method = args.method
    fastq = args.fastq
    fastq_out = args.fastq_out
    cb_count = args.cell_barcode_count
    num_thread = args.num_thread

    scumi.format_fastq(*fastq, config=config, method=method, fastq_out=fastq_out,
                       cb_count=cb_count, num_thread=num_thread)


def tag_bam(args):
    bam = args.bam
    gtf = args.gtf

    annotate_multi_mapping = args.annotate_multi_mapping
    annotate_intron = args.annotate_intron
    strand = args.strand
    num_thread = args.num_thread
    featureCounts = args.featureCounts

    scumi.annotate_bam(bam, gtf=gtf, featureCounts=featureCounts,
                       annotate_multi_mapping=annotate_multi_mapping,
                       strand=strand, annotate_intron=annotate_intron,
                       num_thread=num_thread)


def count_umi(args):
    cb = args.cb
    bam = args.bam
    molecular_info_h5 = args.molecular_info_h5
    gtf = args.gtf
    cb_count = args.cell_barcode_count
    feature_tag = args.feature_tag
    expect_cell = args.expect_cell
    force_cell = args.force_cell
    depth_threshold = args.depth_threshold
    cell_barcode_whitelist = args.cell_barcode_whitelist

    scumi.count_feature(*cb, bam=bam, molecular_info_h5=molecular_info_h5, gtf=gtf,
                        cb_count=cb_count, feature_tag=feature_tag, expect_cell=expect_cell,
                        force_cell=force_cell, depth_threshold=depth_threshold,
                        cell_barcode_whitelist=cell_barcode_whitelist)


def down_sample(args):
    molecular_info = args.molecular_info
    total_read = args.total_read

    total_cell = args.total_cell
    mean_read = args.mean_read
    out_prefix = args.out_prefix
    depth_threshold = args.depth_threshold
    seed = args.seed

    scumi.down_sample(molecular_info, total_read=total_read,
                      total_cell=total_cell,
                      mean_read=mean_read, out_prefix=out_prefix,
                      depth_threshold=depth_threshold,
                      seed=seed)
