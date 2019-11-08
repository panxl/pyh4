"""
Wrapper for the H4 C code
"""
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext as _build_ext


class build_ext(_build_ext):
    def get_ext_filename(self, ext_name):
        return ext_name + '.so'


short_description = __doc__.split("\n")

try:
    with open("README.md", "r") as fh:
        long_description = fh.read()
except:
    long_description = "\n".join(short_description[2:])

setup(
    name="pyh4",
    version="1.0",
    author="Xiaoliang Pan",
    author_email="panxl@panxl.net",
    description=short_description[0],
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/panxl/pyh4",
    packages=["pyh4"],
    ext_modules=[Extension("pyh4/libh4", ["src/h_bonds4.c"])],
    cmdclass={'build_ext': build_ext},
    zip_safe=False,
)
