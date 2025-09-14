"""
Phylogenetic Tree Analyzer

A comprehensive Python package for analyzing evolutionary relationships among species.
Supports sequence alignment, phylogenetic tree construction, and visualization.
"""

__version__ = "1.0.0"
__author__ = "Phylogenetic Tree Analyzer Team"
__email__ = "contact@phyloanalyzer.com"

from .sequence_parser import SequenceParser
from .alignment_manager import AlignmentManager
from .tree_builder import TreeBuilder
from .tree_visualizer import TreeVisualizer
from .main import main

__all__ = [
    'SequenceParser',
    'AlignmentManager', 
    'TreeBuilder',
    'TreeVisualizer',
    'main'
]
