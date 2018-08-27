import gzip
import itertools
from io import BufferedReader, TextIOWrapper


def read_fastq(fastq):
    """
    Return fastq records
    """
    if fastq.endswith('gz'):
        fastq_file = TextIOWrapper(BufferedReader(gzip.open(fastq, mode='rb')))
    else:
        fastq_file = open(fastq, mode='rt')

    while True:
        element = ''.join(itertools.islice(fastq_file, 4))
        if element is not '':
            yield element
        else:
            break
    fastq_file.close()

    return element


def write_fastq(fastq):
    """
    Open a fastq file for writing
    """
    fastq_file = None
    if fastq:
        if fastq.endswith('gz'):
            fastq_file = gzip.open(fastq, mode='wb')
        else:
            fastq_file = open(fastq, mode='wt')

    return fastq_file
