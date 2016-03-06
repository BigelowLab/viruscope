from setuptools import setup


_locals = {}
with open('viralscan.py') as fp:
    exec(fp.read(), None, _locals)
version = _locals['__version__']

setup(
    name='viralscan',
    version=version,
    url='https://github.com/BigelowLab/viralscan',
    license='',
    author='Joe Brown', 'Ben Tupper'
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
