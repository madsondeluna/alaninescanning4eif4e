from setuptools import setup, find_packages

setup(
    name="rosetta-alanine-scanning",
    version="0.1.0",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    install_requires=[
        "click>=8.1.7",
        "biopython>=1.81",
        "pandas>=2.0.0",
        "numpy>=1.24.0",
        "matplotlib>=3.7.0",
        "seaborn>=0.12.0",
        "pyyaml>=6.0",
        "rich>=13.0.0",
    ],
    entry_points={
        "console_scripts": [
            "rosetta-scan=rosetta_scan.cli:main",
        ],
    },
    python_requires=">=3.8",
    author="Madson Luna",
    description="Elegant Rosetta Flex ddG protocol for alanine scanning",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
