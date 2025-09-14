"""
Setup script for Phylogenetic Tree Analyzer package.
"""

from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

with open("requirements.txt", "r", encoding="utf-8") as fh:
    requirements = [line.strip() for line in fh if line.strip() and not line.startswith("#")]

setup(
    name="phylogenetic-tree-analyzer",
    version="1.0.0",
    author="Phylogenetic Tree Analyzer Team",
    author_email="contact@phyloanalyzer.com",
    description="A comprehensive Python package for analyzing evolutionary relationships among species",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/phylogenetic-tree-analyzer",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Information Analysis",
    ],
    python_requires=">=3.8",
    install_requires=requirements,
    entry_points={
        "console_scripts": [
            "phylo-analyzer=phylogenetic_tree_analyzer.main:main",
        ],
    },
    include_package_data=True,
    package_data={
        "phylogenetic_tree_analyzer": ["*.txt", "*.md"],
    },
    keywords="phylogenetics, bioinformatics, evolution, tree, alignment, sequences",
    project_urls={
        "Bug Reports": "https://github.com/yourusername/phylogenetic-tree-analyzer/issues",
        "Source": "https://github.com/yourusername/phylogenetic-tree-analyzer",
        "Documentation": "https://github.com/yourusername/phylogenetic-tree-analyzer/blob/main/README.md",
    },
)
