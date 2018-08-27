import sys
import numpy as np
import logging
import time
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s: %(levelname)s: %(message)s')
logger = logging.getLogger(__name__)


def get_cell_whitelist(cb_count,
                       expect_cell=False,
                       force_cell=False,
                       plot_prefix='.'):

    time_whitelist = time.time()
    cell_whitelist = detect_cell_barcode(
        cb_count, expect_cell, force_cell, plot_prefix)
    time_whitelist = time.time() - time_whitelist
    logger.info(f'Detecting whitelist done, taking {time_whitelist/60.0:.3f} minutes')

    # Quality estimation
    if cell_whitelist is None:
        logger.error('No cell barcode detected, aborting!')
        sys.exit()

    return cell_whitelist


# TODO: cell-barcode class
def detect_cell_barcode(cb_count,
                        expect_cell=False,
                        force_cell=False,
                        plot_prefix='.'):

    count = sorted(cb_count.values(), reverse=True)

    if force_cell:
        cell_number = int(expect_cell)
        cell_number = min(cell_number, len(count))
        threshold = count[cell_number]
    else:
        threshold = get_default_cell_barcode(count=count, expect_cell=expect_cell)

    final_barcodes = set([
        x for x, y in cb_count.items() if y > threshold])
    cell_number = len(final_barcodes)

    plot_barcode_detection(count, cell_number, plot_prefix)

    return final_barcodes


def get_default_cell_barcode(count, expect_cell=False):
    log_cumsum_count = np.log10(np.cumsum(count))

    if not expect_cell:
        total_umi = log_cumsum_count[-1] * 0.8
        expect_cell = np.sum([log_cumsum_count <= total_umi])
    else:
        expect_cell = int(expect_cell)

    if expect_cell <= 500:
        q75, q25 = np.percentile(count[:expect_cell], [75, 25])
        umi_threshold = q75 + (q75 - q25) * 1.5
        umi_threshold /= 10.0
    else:
        umi_threshold = np.percentile(count[:expect_cell], 99) / 10.0

    return umi_threshold


def plot_barcode_detection(count, cell_number, plot_prefix):
    # Knee plot
    fig = plt.figure()
    fig1 = fig.add_subplot(1, 1, 1)

    log_cumsum_count = np.log10(np.cumsum(count))
    if cell_number:
        fig1.scatter(x=range(1, cell_number + 1), y=log_cumsum_count[:cell_number],
                     c='black', s=2.5, label='Cells (' + str(cell_number) + ')')
        fig1.scatter(range(cell_number + 1, len(count) + 1), log_cumsum_count[cell_number:],
                     c='silver', s=2.0, label='Background')
    else:
        fig1.scatter(range(1, len(count) + 1), log_cumsum_count,
                     c='silver', s=2.0, label='Background')

    fig1.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

    fig1.set_xlabel('Cell barcode index')
    fig1.set_ylabel('Cumulative counts (log10)')

    fig.savefig(f'{plot_prefix}_cell_barcode_knee_plot.png',
                bbox_inches='tight', dpi=600)

    # Count plot
    fig = plt.figure()
    fig1 = fig.add_subplot(1, 1, 1)

    if cell_number:
        fig1.scatter(x=range(1, len(count) + 1), y=count,
                     c='white', s=2.0)
        fig1.scatter(x=range(1, cell_number + 1), y=count[:cell_number],
                     c='black', s=2.5, label='Cells (' + str(cell_number) + ')')
        fig1.scatter(x=range(cell_number + 1, len(count) + 1), y=count[cell_number:],
                     c='silver', s=2.0, label='Background')
        fig1.axvline(x=cell_number, linewidth=1.0, c='gray')
    else:
        fig1.scatter(x=range(1, len(count) + 1), y=count,
                     c='silver', s=2.0)

    fig1.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

    fig1.loglog()
    fig1.set_xlim(0.95, len(count) * 1.25)
    fig1.set_xlabel('Cell barcode index')
    fig1.set_ylabel('Read count')

    fig.savefig(f'{plot_prefix}_cell_barcode_count_plot.png',
                bbox_inches='tight', dpi=600)
