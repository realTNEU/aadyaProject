"""
Tree Visualizer Module

This module handles visualization of phylogenetic trees using matplotlib
and other plotting libraries.
"""

import logging
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from typing import List, Dict, Optional, Tuple, Union
import numpy as np
from Bio import Phylo
from Bio.Phylo.BaseTree import Tree, Clade
import seaborn as sns
from matplotlib.patches import FancyBboxPatch
import warnings

# Suppress matplotlib warnings
warnings.filterwarnings('ignore', category=UserWarning, module='matplotlib')


class TreeVisualizer:
    """
    A class for visualizing phylogenetic trees with various customization options.
    
    Supports different tree layouts, color schemes, and annotation options.
    """
    
    def __init__(self, logger: Optional[logging.Logger] = None):
        """
        Initialize the TreeVisualizer.
        
        Args:
            logger: Optional logger instance. If None, creates a default logger.
        """
        self.logger = logger or logging.getLogger(__name__)
        self.default_figsize = (12, 8)
        self.default_dpi = 300
        
        # Set up matplotlib style
        plt.style.use('default')
        sns.set_palette("husl")
        
    def plot_tree(self, tree: Phylo.BaseTree.Tree,
                  title: str = "Phylogenetic Tree",
                  figsize: Tuple[int, int] = None,
                  dpi: int = None,
                  show_branch_lengths: bool = True,
                  show_bootstrap_values: bool = True,
                  color_scheme: str = 'default',
                  node_size: int = 300,
                  font_size: int = 10,
                  output_path: Optional[str] = None) -> plt.Figure:
        """
        Plot a phylogenetic tree with various customization options.
        
        Args:
            tree: Phylogenetic tree object to plot
            title: Title for the plot
            figsize: Figure size tuple (width, height)
            dpi: Dots per inch for the plot
            show_branch_lengths: Whether to show branch lengths
            show_bootstrap_values: Whether to show bootstrap values
            color_scheme: Color scheme for the tree ('default', 'viridis', 'plasma', 'custom')
            node_size: Size of internal nodes
            font_size: Font size for labels
            output_path: Path to save the plot (optional)
            
        Returns:
            Matplotlib figure object
        """
        try:
            figsize = figsize or self.default_figsize
            dpi = dpi or self.default_dpi
            
            # Create figure and axis
            fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
            
            # Plot the tree
            self._draw_tree(ax, tree, show_branch_lengths, show_bootstrap_values,
                          color_scheme, node_size, font_size)
            
            # Customize the plot
            ax.set_title(title, fontsize=16, fontweight='bold', pad=20)
            ax.set_xlabel('Branch Length', fontsize=12)
            ax.set_ylabel('Taxa', fontsize=12)
            
            # Remove spines and ticks for cleaner look
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['left'].set_visible(False)
            ax.tick_params(left=False, labelleft=False)
            
            # Adjust layout
            plt.tight_layout()
            
            if output_path:
                self.save_plot(fig, output_path)
                self.logger.info(f"Tree plot saved to {output_path}")
            
            self.logger.info("Tree visualization completed")
            return fig
            
        except Exception as e:
            self.logger.error(f"Error plotting tree: {str(e)}")
            raise
    
    def _draw_tree(self, ax: plt.Axes, tree: Phylo.BaseTree.Tree,
                   show_branch_lengths: bool, show_bootstrap_values: bool,
                   color_scheme: str, node_size: int, font_size: int) -> None:
        """
        Draw the phylogenetic tree on the given axis.
        
        Args:
            ax: Matplotlib axis to draw on
            tree: Phylogenetic tree object
            show_branch_lengths: Whether to show branch lengths
            show_bootstrap_values: Whether to show bootstrap values
            color_scheme: Color scheme for the tree
            node_size: Size of internal nodes
            font_size: Font size for labels
        """
        try:
            # Get tree coordinates
            coords = self._get_tree_coordinates(tree)
            
            # Set up color palette
            colors = self._get_color_palette(color_scheme, len(coords))
            
            # Draw branches
            self._draw_branches(ax, tree, coords, colors, show_branch_lengths, font_size)
            
            # Draw nodes
            self._draw_nodes(ax, tree, coords, colors, node_size, show_bootstrap_values, font_size)
            
            # Draw labels
            self._draw_labels(ax, tree, coords, font_size)
            
        except Exception as e:
            self.logger.error(f"Error drawing tree: {str(e)}")
            raise
    
    def _get_tree_coordinates(self, tree: Phylo.BaseTree.Tree) -> Dict[Clade, Tuple[float, float]]:
        """
        Calculate coordinates for tree nodes.
        
        Args:
            tree: Phylogenetic tree object
            
        Returns:
            Dictionary mapping clades to (x, y) coordinates
        """
        try:
            # This is a simplified coordinate calculation
            # In practice, you might want to use more sophisticated layout algorithms
            
            coords = {}
            terminals = tree.get_terminals()
            internals = tree.get_nonterminals()
            
            # Calculate y-coordinates for terminals
            terminal_y = {}
            for i, terminal in enumerate(terminals):
                terminal_y[terminal] = i
            
            # Calculate x-coordinates (branch lengths from root)
            def calculate_x_coords(clade, parent_x=0.0):
                if clade.branch_length is not None:
                    x = parent_x + clade.branch_length
                else:
                    x = parent_x + 0.1  # Default branch length
                
                coords[clade] = (x, 0)  # Will update y-coordinates later
                
                for child in clade.clades:
                    calculate_x_coords(child, x)
            
            # Start from root
            if tree.root:
                calculate_x_coords(tree.root)
            
            # Update y-coordinates for internal nodes
            def update_y_coords(clade):
                if clade.is_terminal():
                    coords[clade] = (coords[clade][0], terminal_y[clade])
                else:
                    # Average y-coordinate of children
                    child_y_coords = [coords[child][1] for child in clade.clades]
                    avg_y = np.mean(child_y_coords)
                    coords[clade] = (coords[clade][0], avg_y)
                    
                    for child in clade.clades:
                        update_y_coords(child)
            
            if tree.root:
                update_y_coords(tree.root)
            
            return coords
            
        except Exception as e:
            self.logger.error(f"Error calculating tree coordinates: {str(e)}")
            return {}
    
    def _get_color_palette(self, color_scheme: str, num_colors: int) -> List[str]:
        """
        Get color palette for the tree.
        
        Args:
            color_scheme: Name of the color scheme
            num_colors: Number of colors needed
            
        Returns:
            List of color strings
        """
        if color_scheme == 'viridis':
            return plt.cm.viridis(np.linspace(0, 1, num_colors))
        elif color_scheme == 'plasma':
            return plt.cm.plasma(np.linspace(0, 1, num_colors))
        elif color_scheme == 'custom':
            return ['#FF6B6B', '#4ECDC4', '#45B7D1', '#96CEB4', '#FECA57', 
                   '#FF9FF3', '#54A0FF', '#5F27CD', '#00D2D3', '#FF9F43']
        else:  # default
            return ['#2E86AB', '#A23B72', '#F18F01', '#C73E1D', '#7209B7',
                   '#0B6623', '#FF6F00', '#E91E63', '#3F51B5', '#009688']
    
    def _draw_branches(self, ax: plt.Axes, tree: Phylo.BaseTree.Tree,
                      coords: Dict[Clade, Tuple[float, float]], colors: List[str],
                      show_branch_lengths: bool, font_size: int) -> None:
        """
        Draw tree branches.
        
        Args:
            ax: Matplotlib axis to draw on
            tree: Phylogenetic tree object
            coords: Dictionary of node coordinates
            colors: List of colors
            show_branch_lengths: Whether to show branch lengths
            font_size: Font size for labels
        """
        try:
            def draw_branch(clade, color_index=0):
                if clade.is_terminal():
                    return
                
                x, y = coords[clade]
                color = colors[color_index % len(colors)]
                
                for child in clade.clades:
                    child_x, child_y = coords[child]
                    
                    # Draw horizontal branch
                    ax.plot([x, child_x], [y, y], color=color, linewidth=2, alpha=0.7)
                    
                    # Draw vertical branch
                    ax.plot([child_x, child_x], [y, child_y], color=color, linewidth=2, alpha=0.7)
                    
                    # Add branch length label if requested
                    if show_branch_lengths and child.branch_length is not None:
                        mid_x = (x + child_x) / 2
                        ax.text(mid_x, y - 0.1, f'{child.branch_length:.3f}', 
                               ha='center', va='top', fontsize=font_size-2, alpha=0.8)
                    
                    # Recursively draw child branches
                    draw_branch(child, color_index + 1)
            
            if tree.root:
                draw_branch(tree.root)
                
        except Exception as e:
            self.logger.error(f"Error drawing branches: {str(e)}")
    
    def _draw_nodes(self, ax: plt.Axes, tree: Phylo.BaseTree.Tree,
                   coords: Dict[Clade, Tuple[float, float]], colors: List[str],
                   node_size: int, show_bootstrap_values: bool, font_size: int) -> None:
        """
        Draw tree nodes.
        
        Args:
            ax: Matplotlib axis to draw on
            tree: Phylogenetic tree object
            coords: Dictionary of node coordinates
            colors: List of colors
            node_size: Size of internal nodes
            show_bootstrap_values: Whether to show bootstrap values
            font_size: Font size for labels
        """
        try:
            for i, (clade, (x, y)) in enumerate(coords.items()):
                if clade.is_terminal():
                    # Terminal nodes
                    ax.scatter(x, y, s=node_size//2, c=colors[i % len(colors)], 
                             marker='o', edgecolors='black', linewidth=1, zorder=3)
                else:
                    # Internal nodes
                    ax.scatter(x, y, s=node_size, c=colors[i % len(colors)], 
                             marker='s', edgecolors='black', linewidth=1, zorder=3)
                    
                    # Add bootstrap value if available
                    if show_bootstrap_values and hasattr(clade, 'confidence') and clade.confidence is not None:
                        ax.text(x, y + 0.2, f'{clade.confidence:.0f}', 
                               ha='center', va='bottom', fontsize=font_size-2, 
                               fontweight='bold', bbox=dict(boxstyle='round,pad=0.2', 
                               facecolor='white', alpha=0.8))
                
        except Exception as e:
            self.logger.error(f"Error drawing nodes: {str(e)}")
    
    def _draw_labels(self, ax: plt.Axes, tree: Phylo.BaseTree.Tree,
                    coords: Dict[Clade, Tuple[float, float]], font_size: int) -> None:
        """
        Draw tree labels.
        
        Args:
            ax: Matplotlib axis to draw on
            tree: Phylogenetic tree object
            coords: Dictionary of node coordinates
            font_size: Font size for labels
        """
        try:
            for clade, (x, y) in coords.items():
                if clade.is_terminal() and clade.name:
                    ax.text(x + 0.1, y, clade.name, va='center', fontsize=font_size, 
                           fontweight='bold')
                
        except Exception as e:
            self.logger.error(f"Error drawing labels: {str(e)}")
    
    def plot_tree_comparison(self, trees: List[Phylo.BaseTree.Tree],
                           tree_names: List[str],
                           title: str = "Tree Comparison",
                           figsize: Tuple[int, int] = None) -> plt.Figure:
        """
        Plot multiple trees side by side for comparison.
        
        Args:
            trees: List of phylogenetic tree objects
            tree_names: List of names for the trees
            title: Title for the plot
            figsize: Figure size tuple
            
        Returns:
            Matplotlib figure object
        """
        try:
            if len(trees) != len(tree_names):
                raise ValueError("Number of trees must match number of tree names")
            
            figsize = figsize or (16, 6)
            fig, axes = plt.subplots(1, len(trees), figsize=figsize)
            
            if len(trees) == 1:
                axes = [axes]
            
            for i, (tree, name) in enumerate(zip(trees, tree_names)):
                self._draw_tree(axes[i], tree, True, True, 'default', 300, 10)
                axes[i].set_title(name, fontsize=14, fontweight='bold')
            
            plt.suptitle(title, fontsize=16, fontweight='bold')
            plt.tight_layout()
            
            self.logger.info("Tree comparison plot completed")
            return fig
            
        except Exception as e:
            self.logger.error(f"Error plotting tree comparison: {str(e)}")
            raise
    
    def save_plot(self, fig: plt.Figure, output_path: str, 
                  format: str = 'png', dpi: int = None) -> None:
        """
        Save a matplotlib figure to a file.
        
        Args:
            fig: Matplotlib figure object
            output_path: Path where to save the plot
            format: Output format ('png', 'pdf', 'svg', 'eps')
            dpi: Dots per inch for raster formats
        """
        try:
            dpi = dpi or self.default_dpi
            fig.savefig(output_path, format=format, dpi=dpi, bbox_inches='tight')
            self.logger.info(f"Plot saved to {output_path}")
        except Exception as e:
            self.logger.error(f"Error saving plot: {str(e)}")
            raise
    
    def create_summary_plot(self, alignment_stats: Dict, tree_stats: Dict,
                          output_path: Optional[str] = None) -> plt.Figure:
        """
        Create a summary plot showing alignment and tree statistics.
        
        Args:
            alignment_stats: Dictionary of alignment statistics
            tree_stats: Dictionary of tree statistics
            output_path: Path to save the plot (optional)
            
        Returns:
            Matplotlib figure object
        """
        try:
            fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 10))
            
            # Alignment length distribution
            if 'alignment_length' in alignment_stats:
                ax1.bar(['Alignment Length'], [alignment_stats['alignment_length']], 
                       color='skyblue', alpha=0.7)
                ax1.set_title('Alignment Length', fontweight='bold')
                ax1.set_ylabel('Length (bp)')
            
            # Conservation percentage
            if 'conservation_percentage' in alignment_stats:
                ax2.bar(['Conservation'], [alignment_stats['conservation_percentage']], 
                       color='lightgreen', alpha=0.7)
                ax2.set_title('Sequence Conservation', fontweight='bold')
                ax2.set_ylabel('Percentage (%)')
                ax2.set_ylim(0, 100)
            
            # Tree statistics
            if 'num_terminals' in tree_stats:
                ax3.bar(['Terminals'], [tree_stats['num_terminals']], 
                       color='lightcoral', alpha=0.7)
                ax3.set_title('Number of Taxa', fontweight='bold')
                ax3.set_ylabel('Count')
            
            # Tree length
            if 'tree_length' in tree_stats:
                ax4.bar(['Tree Length'], [tree_stats['tree_length']], 
                       color='gold', alpha=0.7)
                ax4.set_title('Total Tree Length', fontweight='bold')
                ax4.set_ylabel('Length')
            
            plt.suptitle('Phylogenetic Analysis Summary', fontsize=16, fontweight='bold')
            plt.tight_layout()
            
            if output_path:
                self.save_plot(fig, output_path)
            
            self.logger.info("Summary plot created")
            return fig
            
        except Exception as e:
            self.logger.error(f"Error creating summary plot: {str(e)}")
            raise
