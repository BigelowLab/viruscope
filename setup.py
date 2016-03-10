import io
from os.path import join, dirname
from setuptools import setup


def get_version(relpath):
  '''Read version info from a file without importing it'''
  for line in io.open(join(dirname(__file__), relpath), encoding='cp437'):
    if '__version__' in line:
      if '"' in line:
        # __version__ = "0.9"
        return line.split('"')[1]
      elif "'" in line:
        return line.split("'")[1]

setup(
    name='viralscan',
    version=get_version("viralscan.py"),
    url='https://github.com/BigelowLab/viralscan',
    license='',
    author='Joe Brown, Ben Tupper',
    author_email='jmbrown@bigelow.org',
    description='',
    long_description=__doc__,
    py_modules=['viralscan'],
    install_requires=[],
    entry_points='''
        [console_scripts]
        viral-scan=viralscan:main
    '''
)
