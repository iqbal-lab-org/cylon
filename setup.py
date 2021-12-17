import glob
from setuptools import setup, find_packages


with open("requirements.txt") as f:
    install_requires = [x.rstrip() for x in f]

setup(
    name="viridian",
    version="0.1.0",
    description="Virus assembler from amplicon sequencing reads",
    packages=find_packages(exclude=["tests"]),
    author="Martin Hunt",
    author_email="mhunt@ebi.ac.uk",
    url="https://github.com/iqbal-lab-org/viridian",
    test_suite="nose.collector",
    tests_require=["pytest"],
    entry_points={"console_scripts": ["viridian = viridian.__main__:main"]},
    install_requires=install_requires,
    license="MIT",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Programming Language :: Python :: 3 :: Only",
        "License :: OSI Approved :: MIT License",
    ],
)
