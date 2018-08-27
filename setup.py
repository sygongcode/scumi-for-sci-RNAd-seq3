from distutils.core import setup

setup(name='scumi',
      version='0.1.0',
      description='Summarizing scRNA-seq data with unique molecular identifier (UMI)',
      author='Jiarui Ding',
      author_email='jiarui.ding@gmail.com',
      package_dir={'': 'lib'},
      packages=['scumi'],
      scripts=['scumi']
      )
