try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

setup(name='scumi',
      version='0.1.0',
      description='Summarizing scRNA-seq data with unique molecular identifier (UMI)',
      author='Jiarui Ding, Sean K. Simmons, Joshua Z. Levin, Aviv Regev',
      author_email='jding@broadinstitute.org',
      package_dir={'': 'lib'},
      packages=['scumi'],
      scripts=['scumi'],
      install_requires=[
          'numpy >= 1.13.1',
          'scipy >= 0.19.1',
          'numba >= 0.37.0',
          'cytoolz >= 0.8.2',
          'pyyaml >= 3.12',
          'regex >= 2.4.127',
          'tables >= 3.4.4',
          'pandas >= 0.21.1',
          'pysam >= 0.11.2.2',
          'matplotlib > 2.0.2'
      ]
      )
