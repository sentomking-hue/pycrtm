import setuptools
from skbuild import setup

requires=['numpy']

setup(
    name="pycrtm",
    version='0.0.1',
    description='Python wrapper for the CRTM.',
    author='Bryan Karpowicz',
    requires=requires,
    packages=['pycrtm'],
    py_modules=['crtm_io', 'pyCRTM']
)
