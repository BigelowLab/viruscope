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
    name='viruscope',
    version=get_version("viruscope.py"),
    url='https://github.com/BigelowLab/viruscope',
    license='',
    author='Joe Brown, Ben Tupper, Julia Brown',
    author_email='jmbrown@bigelow.org',
    description='',
    long_description=__doc__,
    py_modules=['viruscope'],
    packages=['viruscope']
    install_requires=[],
    entry_points='''
        [console_scripts]
        viruscope=viruscope:main
    '''
)
