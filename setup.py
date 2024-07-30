from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="sesi_prep",  
    version="1.0",
    author="Shi Yao",
    author_email="yaoshi229@gmail.com",
    description="A package for prepare seismic data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/YaoShi-igg/seis_prep",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    install_requires=[
        "numpy>=1.24.3",
        "matplotlib>=3.7.2",
        "obspy>=1.4.0"
    ],
)

