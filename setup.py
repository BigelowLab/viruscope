from setuptools import setup


setup(
    name='viralscan',
    version='0.1.9',
    url='https://github.com/BigelowLab/viralscan',
    license='',
    author='Joe Brown',
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
