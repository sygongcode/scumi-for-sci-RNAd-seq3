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
