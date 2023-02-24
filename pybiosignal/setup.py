from setuptools import setup

from pybiosignal import __version__

setup(
    name='my_pip_package',
    version=__version__,

    url='https://github.com/donnubio/gitlib',

    author='artembio2',
    author_email='artembio2@gmail.com',

    py_modules=['pybiosignal'],
)
