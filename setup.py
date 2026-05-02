from setuptools import setup, find_packages

setup(
    name="alleleselect",
    version="1.0.0",
    author="Angie Xiu",
    author_email="angie.xiu27@gmail.com",
    description="Allele-Selective ASO Design Pipeline for CACNA1A Gain-of-Function Mutations",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/axshoe/alleleselect",
    packages=find_packages(),
    python_requires=">=3.9",
    install_requires=[
        "requests>=2.28.0",
        "biopython>=1.81",
        "numpy>=1.24.0",
        "scipy>=1.10.0",
    ],
    extras_require={
        "dev": ["pytest>=7.0"],
    },
    entry_points={
        "console_scripts": [
            "alleleselect=alleleselect.cli:main",
        ],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)
