# Phylogenetic Tree Analyzer

A comprehensive Python package for analyzing evolutionary relationships among species through phylogenetic tree construction and visualization.

## Features

- **Sequence Parsing**: Support for FASTA and GenBank file formats
- **Multiple Sequence Alignment**: MUSCLE, ClustalW, and simple alignment methods
- **Phylogenetic Tree Construction**: Neighbor-Joining, UPGMA, and Maximum Parsimony methods
- **Tree Visualization**: High-quality tree plots with customizable styling
- **Bootstrap Analysis**: Support for bootstrap replicates to assess tree reliability
- **Comprehensive Logging**: Detailed logging and error handling
- **Command-line Interface**: Easy-to-use CLI with extensive options
- **Modular Design**: Well-structured, documented code for easy extension

## Installation

### Prerequisites

- Python 3.8 or higher
- pip package manager

### Install Dependencies

```bash
pip install -r requirements.txt
```

### Install the Package

```bash
pip install -e .
```

## Quick Start

### Basic Usage

```bash
# Analyze a FASTA file with default parameters
python -m phylogenetic_tree_analyzer data/your_sequences.fasta -o output/

# Analyze with custom parameters
python -m phylogenetic_tree_analyzer data/your_sequences.fasta -o output/ \
    --alignment-method muscle \
    --tree-method neighbor_joining \
    --bootstrap-replicates 1000 \
    --log-level DEBUG
```

### Example Workflow

1. **Prepare your sequences**: Place your FASTA file in the `data/` folder
2. **Run analysis**: Execute the command with your input file
3. **View results**: Check the `output/` folder for:
   - Aligned sequences (`alignments/`)
   - Phylogenetic tree (`trees/`)
   - Tree visualization (`plots/`)
   - Analysis report (`analysis_report.txt`)
   - Log file (`logs/`)

## Command-line Options

### Required Arguments
- `input_file`: Path to input sequence file (FASTA or GenBank format)
- `-o, --output-dir`: Output directory for results

### Analysis Parameters
- `--alignment-method`: Choose from `muscle`, `clustalw`, or `simple` (default: muscle)
- `--tree-method`: Choose from `neighbor_joining`, `upgma`, or `maximum_parsimony` (default: neighbor_joining)
- `--substitution-model`: Choose from `identity`, `blastn`, or `trans` (default: identity)
- `--bootstrap-replicates`: Number of bootstrap replicates (default: 100)
- `--min-sequence-length`: Minimum sequence length filter (default: 50)

### Output Options
- `--log-level`: Choose from `DEBUG`, `INFO`, `WARNING`, `ERROR` (default: INFO)

## Project Structure

```
phylogenetic_tree_analyzer/
├── phylogenetic_tree_analyzer/
│   ├── __init__.py
│   ├── main.py                 # Command-line interface
│   ├── sequence_parser.py      # Sequence parsing module
│   ├── alignment_manager.py    # Alignment algorithms
│   ├── tree_builder.py         # Phylogenetic tree construction
│   └── tree_visualizer.py      # Tree visualization
├── data/                       # Input sequence files
├── output/                     # Analysis results
│   ├── alignments/            # Aligned sequences
│   ├── trees/                 # Phylogenetic trees
│   ├── plots/                 # Tree visualizations
│   └── logs/                  # Log files
├── requirements.txt           # Dependencies
└── README.md                  # This file
```

## Supported File Formats

### Input Formats
- **FASTA** (`.fasta`, `.fa`): Standard sequence format
- **GenBank** (`.genbank`, `.gb`, `.gbk`): Annotated sequence format

### Output Formats
- **Alignments**: FASTA, Clustal, PHYLIP
- **Trees**: Newick, NEXUS, PhyloXML
- **Plots**: PNG, PDF, SVG, EPS

## API Usage

You can also use the package programmatically:

```python
from phylogenetic_tree_analyzer import (
    SequenceParser, AlignmentManager, TreeBuilder, TreeVisualizer
)

# Parse sequences
parser = SequenceParser()
sequences = parser.parse_file('data/sequences.fasta')

# Align sequences
aligner = AlignmentManager()
alignment = aligner.multiple_sequence_align(sequences, method='muscle')

# Build tree
tree_builder = TreeBuilder()
tree = tree_builder.build_tree(alignment, method='neighbor_joining')

# Visualize tree
visualizer = TreeVisualizer()
fig = visualizer.plot_tree(tree, title="My Phylogenetic Tree")
```

## Algorithm Details

### Alignment Methods
- **MUSCLE**: Fast and accurate multiple sequence alignment
- **ClustalW**: Classic progressive alignment algorithm
- **Simple**: Basic pairwise alignment for small datasets

### Tree Construction Methods
- **Neighbor-Joining**: Fast distance-based method, good for large datasets
- **UPGMA**: Unweighted Pair Group Method with Arithmetic Mean
- **Maximum Parsimony**: Character-based method (simplified implementation)

### Substitution Models
- **Identity**: Simple identity matrix
- **BLASTN**: Nucleotide substitution matrix
- **TRANS**: Transition/transversion matrix

## Troubleshooting

### Common Issues

1. **MUSCLE not found**: Ensure MUSCLE is installed and in PATH
2. **ClustalW not found**: Install ClustalW2 and ensure it's in PATH
3. **Memory issues**: Reduce bootstrap replicates or use simpler methods
4. **File format errors**: Ensure input files are valid FASTA/GenBank format

### Getting Help

- Check the log files in `output/logs/` for detailed error messages
- Use `--log-level DEBUG` for verbose output
- Ensure all dependencies are properly installed

## Contributing

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Add tests if applicable
5. Submit a pull request

## License

This project is licensed under the MIT License.

## Citation

If you use this tool in your research, please cite:

```
Phylogenetic Tree Analyzer v1.0.0
A comprehensive tool for phylogenetic analysis
https://github.com/yourusername/phylogenetic-tree-analyzer
```

## Acknowledgments

- Built with [Biopython](https://biopython.org/)
- Visualization powered by [Matplotlib](https://matplotlib.org/)
- Inspired by phylogenetic analysis best practices
