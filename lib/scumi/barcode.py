import sys
import itertools
import numpy as np
import logging
from collections import defaultdict
from .util import compute_edit_distance

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s: %(levelname)s: %(message)s')
logger = logging.getLogger(__name__)


class MergeHash(object):
    def __init__(self, barcode1, barcode2, rev=False):
        merge_hash = defaultdict(set)
        if rev is True:
            for idx in range(len(barcode1)):
                dt_cb = barcode1[idx]
                merge_hash[dt_cb].update([dt_cb])
                random_cb = barcode2[idx]
                merge_hash[dt_cb].update([random_cb])
        else:
            for idx in range(len(barcode1)):
                dt_cb = barcode1[idx]
                merge_hash[dt_cb].update([dt_cb])
                random_cb = barcode2[idx]
                merge_hash[random_cb].update([dt_cb])
        self.hash = merge_hash

    def __getitem__(self, barcode):
        corrected_barcode = self.hash[barcode]
        if len(corrected_barcode) != 1:
            return None
        else:
            return list(corrected_barcode)[0]


# TODO barcode quality check
class ErrorBarcodeHash(object):
    def __init__(self, barcodes, edit_distance=1):
        if not isinstance(edit_distance, int) or edit_distance < 0:
            logger.error('edit_distance should be an integer and non-negative, '
                         'aborting.')
            sys.exit(-1)

        if isinstance(barcodes, str):
            barcodes = [barcodes]
        self.barcode_len = len(barcodes[0])

        if len(barcodes) <= 0 or self.barcode_len <= 0:
            logger.error('The number of barcodes and the barcode length '
                         'should be greater than 1')

        self.hash = defaultdict(set)
        self.generate_error_hash_table(barcodes, edit_distance)

    def __getitem__(self, barcode):
        corrected_barcode = self.hash[barcode]
        if len(corrected_barcode) != 1:
            return None
        else:
            return list(corrected_barcode)[0]

    def generate_error_hash_table(self, barcodes, edit_distance=1):
        if self.barcode_len <= edit_distance:
            logger.error('edit_distance should be less than barcode length, '
                         'aborting.')
            sys.exit(-1)

        index, value = \
            ErrorBarcodeHash.generate_index_value(self.barcode_len, edit_distance)

        collision_barcode = 0
        for barcode in barcodes:
            self.hash[barcode].update([barcode])

            for x in ErrorBarcodeHash.generate_error_barcode(barcode, index, value):
                if x not in barcodes:
                    self.hash[x].update([barcode])
                elif x != barcode:
                    collision_barcode += 1
                    logger.info(f'Warning: barcode {x} and {barcode} are '
                                f'within edit distance of {edit_distance}')

        logger.info(f'The number of potential collision barcodes '
                    f'(edit distance {edit_distance}): {collision_barcode}\n')

    @staticmethod
    def generate_index_value(barcode_len, edit_distance=1):
        alphabet = 'ACGTN'

        index = itertools.combinations(range(barcode_len), edit_distance)
        value = itertools.product(alphabet, repeat=edit_distance)

        return list(index), list(value)

    @staticmethod
    def generate_error_barcode(barcode, index, value):

        def mutate_barcode(idx, v):
            bc = list(barcode)
            for i, i_v in enumerate(idx):
                bc[i_v] = v[i]
            return ''.join(bc)

        return {mutate_barcode(idx, v) for idx in index for v in value}


class ErrorBarcodeHashConstraint(ErrorBarcodeHash):
    def __init__(self, barcodes, barcode_hash, edit_distance=1):
        super().__init__(barcodes, edit_distance)

        self.barcode_hash = barcode_hash

    def __getitem__(self, barcode):
        barcode, barcode_full = barcode

        corrected_barcode = list(self.hash[barcode])

        num_candidate = len(corrected_barcode)
        if num_candidate == 1:
            return corrected_barcode[0]
        elif num_candidate == 0:
            return None
        else:
            dist_all = np.zeros(num_candidate, dtype=np.int32)
            for idx, x in enumerate(corrected_barcode):
                barcode_list = list(self.barcode_hash[x])
                edit_dist = [compute_edit_distance(barcode_full.encode('utf-8'),
                                                   x.encode('utf-8'))
                             for x in barcode_list]
                dist_all[idx] = min(edit_dist)

            idx = np.where(dist_all == min(dist_all))[0]
            if len(idx) == 1:
                return corrected_barcode[idx[0]]
            else:
                return None
