---
output:
  html_document: default
  pdf_document: default
---
![scumi](doc/scumi.png)


Summarizing single-cell RNA-sequencing data with unified molecular identifiers
====================

scumi is a flexible Python package to process fastq files generated from different single-cell RNA-sequencing (scRNA-seq) protocols to produce a gene-cell sparse expression matrix for downstream analyses, e.g., discovering cell types and inferring cell lineages. 


## Installation 
scumi has been tested on Python 3.6 on both Linux (Linux version 2.6.32-696.30.1.el6.x86_64) and macOS Sierra (Version 10.12.6). 

Dependencies

* Python 3.6 or above 
* [The STAR aligner](https://github.com/alexdobin/STAR)
* The STAR index file. We use one from [10x Genomics](https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest)
* A gene transfer format (GTF) file. We use one from [10x Genomics](https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest)
* The [Subread package](http://subread.sourceforge.net/), version 1.6.2 or above

scumi can be installed using `python setup.py install`. 
We will simplify the installation using conda packaging.


## Running scumi

Here we used example data from sci-RNA-seq to show how to run scumi. 
The same pipeline can be used to analyze data from other protocols such 10x Chromium, Drop-seq, Seq-Well, CEL-Seq2, InDrop, and SPLiT-seq. 

```bash
# pre-installed software
star_dir=/ahg/regevdata/projects/sc_compare/software/STAR-2.6.1a/bin/Linux_x86_64/
feature_count_dir=/ahg/regevdata/projects/sc_compare/software/subread-1.6.2-Linux-x86_64/bin/

index_dir=/seq/regev_genome_portal/SOFTWARE/10X/refdata-cellranger-1.2.0/refdata-cellranger-GRCh38-1.2.0
star_index=$index_dir/star
gtf=$index_dir/genes/genes.gtf

# The input fastq data
fastq_dir=/ahg/regevdata/projects/sc_compare/hiseq-pbmc/SN0142147/fastq/rename
fastq1=$fastq_dir/CC7W2ANXX.120717_SciSeq-p5_H6.unmapped.1.fastq.gz
fastq2=$fastq_dir/CC7W2ANXX.120717_SciSeq-p5_H6.unmapped.2.fastq.gz

# RT-barcode, the candidate cell barcodes
sciRNAseq_RT_barcode=/ahg/regevdata/projects/sc_compare/doc/RT_384_080717.tsv

# Some important intermediate output files
fastq_out=CC7W2ANXX.120717_SciSeq-p5_H6.unmapped.1.fastq.format_fastq.fastq.gz
cell_barcode_count=CC7W2ANXX.120717_SciSeq-p5_H6.unmapped.1.fastq.format_fastq.cb_count

bam=CC7W2ANXX.120717_SciSeq-p5_H6.unmapped.1.fastq.format_fastq.fastq.star_align.bam
molecular_info_h5=$bam".featureCounts.count_feature.h5"


# Step one: Extract cell barcodes and unified molecular identifiers (UMIs) and 
# put them to the header of the cDNA read
scumi merge_fastq $fastq1 $fastq2 \
  --config /ahg/regevdata/projects/sc_compare/doc/config.yaml  \
  --method sci-RNA-seq-polyT  \
  --fastq_out  $fastq_out \
  --num_thread 4  \
  --cell_barcode_count $cell_barcode_count 
  

# Step two: Align the merged fastq file to the corresponding genome using STAR
$star_dir/STAR \
  --genomeDir $star_index \
  --chimOutType WithinBAM \
  --outSAMunmapped Within \
  --outSAMtype BAM Unsorted \
  --outStd BAM_Unsorted \
  --runThreadN 4 \
  --readFilesCommand zcat  \
  --readFilesIn $fastq_out \
  --outFileNamePrefix $fastq_out \
  --outFilterMultimapNmax 20  \
  --outFilterMismatchNoverLmax 0.06  \
  --limitOutSJcollapsed 2000000  \
  --twopassMode Basic  \
  --limitIObufferSize 400000000 > $bam


# Step three: Using featureCounts to annotate each alignment with a gene tag, XT:Z
scumi tag_bam $bam \
  --gtf $gtf  \
  --featureCounts $feature_count_dir"/featureCounts" \
  --annotate_multi_mapping 


# Step four: Counting the number of UMIs per gene per cell 
scumi count_umi \
  --bam $bam".featureCounts.bam" \
  --molecular_info_h5  $molecular_info_h5 \
  --cell_barcode_count $cell_barcode_count  \
  --gtf $gtf  \
  --expect_cell 40 \
  $sciRNAseq_RT_barcode 
```

### Specifying the read structures in the configration file (step one)
After demultiplexing, the sci-RNA-seq data have two reads.
Read1 has 8bp UMI sequences followed by 10bp cell barcodes. 
For the test data, there are poly-T sequences after the cell barcodes in read1. 
scumi will look at the poly-T sequences from base 19 to base 23 (as specified in the following configration file), and will discard the reads with more than one non-T base. 
If you do not want to filter reads based on the poly-T sequences, you can set the parameter `--method sci-RNA-seq` because, 
as can see from the below configuration file `/ahg/regevdata/projects/sc_compare/doc/config.yaml`, there are no poly-T entries in `sci-RNA-seq`. 
Read2 consists of the actual cDNA sequences. 
```bash
sci-RNA-seq: {
    num_read: 2,

    read1: {
      UMI: [1, 8],
      CB1: [9, 18]
    },

    read2: {
      cDNA: []
    }
}

sci-RNA-seq-polyT: {
    num_read: 2,

    read1: {
      UMI: [1, 8],
      CB1: [9, 18],
      polyT: [19, 23]
    },

    read2: {
      cDNA: []
    }
}
```

### Mapping the merged reads using STAR (step two)
STAR has rich parameter setting and the results are typically robust to the parameters. 
We set parameters based on the ENCODE project, e.g., `--outFilterMultimapNmax 20`, and make some adjustments, e.g.,  `--outFilterMismatchNoverLmax 0.06`. 


### Annotate each alignment with a gene tag (step three)
scumi uses featureCounts to tag each alignment. 
To tag multi-mapped reads, set `--annotate_multi_mapping`, 
and to tag intron mapping reads, set `--annotate_intron`. 


### Counting the number of UMIs in each cell (step four)
The annotated bam file by featureCounts can be used as an input for counting the final gene by cell UMI sparse count matrix.
The `--cell_barcode_count` file is generated in step one, which has the number of reads for each cell barcode and the T-frequency information for correcting the bead synthesization errors in Drop-seq data. 
The `--depth_threshold` parameter is used to filter UMIs that have less than a given number of read support.
This feature is useful for some noisy data, e.g., the CEL-Seq2 data that could have barcode switch issues.
The `--expect_cell` parameter specifies the expected number of cells in the input bam file,
and the `--force_cell` parameter forces scumi to ouput the number of cells specified by `--expect_cell`. 
If the set of cell barcode whitelists is known, scumi can only consider these cell barcodes by specifying `--cell_barcode_whitelist=$cell_barcode_whitelist`, 
where `$cell_barcode_whitelist` is an one column text file with the whitelist of cell barcodes (without header).  


For many scRNA-seq protocols (e.g., sci-RNA-seq, 10x chromium, InDrop, and CEL-Seq2), the cell barcodes are drawn from a pool of known candidate cell barcodes.  For example, for the sci-RNA-seq data, the RT-barcodes (cell barcode one) are drawn from 384 candidate cell barcodes. 
Therefore, we can map the observed cell barcodes to the candidate cell barcodes by putting the candidate cell barcodes in a file (`$sciRNAseq_RT_barcode` in our case) and use it as an input to `count_umi`. 
For other protocols such as InDrop, as there are three cell barcodes for each fragment (two cell barcodes plus a sample barcode), 
we can put the corresponding cell barcodes in three files (without headers) and use these three files as inputs to `count_umi`. 

If for some reasons, the candidate cell barcodes for one cell barcode are unknown, you can use `None` as input. For example, assuming that `CB2` is unknown, we can use `$cell_barcode_file_one None $cell_barcode_three` as inputs, where `$cell_barcode_file_one` and `$cell_barcode_file_three` are the candidate cell barcode files for cell barcode one and cell barcode three, respectively.  

```bash
InDrop-polyT: {
  num_read: 4,

  read1: {
    CB1: [1, 8]
  },

  read2: {
    CB2: [1, 8],
    CB2_value: ["CATAACTG", "GGAGGTAA"]
  },

  read3: {
    cDNA: []
  },

  read4: {
    CB3: [1, 8],
    UMI: [9, 14],
    polyT: [15, 19]
  }
}
```
scumi has another useful feature to only extract the reads with given cell barcodes (with edit distances up to one). 
For example, if we know cell barcode two can be either "CATAACTG" or "GGAGGTAA", we can set the value as in the above configuration file. 

# Pipeline management



##

We currently use [bpipe](https://github.com/ssadedin/bpipe) for pipeline management. 
Using bpipe, we can easily build pipelines by connecting different modules (e.g., different steps of the scumi package) to form a pipeline.
Moreover, when some jobs failed, there is no need to re-run the whole pipeline, but started at the steps where the jobs failed. 



