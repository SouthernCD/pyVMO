# coding utf8
import setuptools
from pyvmo.versions import get_versions

with open('README.md') as f:
    LONG_DESCRIPTION = f.read()

setuptools.setup(
    name="pyVMO",
    version=get_versions(),
    author="Yuxing Xu",
    author_email="xuyuxing@mail.kib.ac.cn",
    description="A test python toolkit for variant site analysis",
    long_description=LONG_DESCRIPTION,
    long_description_content_type='text/markdown',
    url="https://github.com/SouthernCD/pyVMO",
    include_package_data = True,

    # entry_points={
    #     "console_scripts": ["PyVMO = pyvmo.cli:main"]
    # },    

    packages=setuptools.find_packages(),

    install_requires=[
        "toolbiox>=0.0.44",
        "scikit-allel>=1.3.7",
        "cyvcf2>=0.30.28",
        "h5py>=3.10.0",
    ],

    python_requires='>=3.5',
)