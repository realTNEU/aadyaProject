"""
Alignment Manager Module

This module handles sequence alignment using various algorithms
provided by Biopython.
"""

import logging
from typing import List, Dict, Optional, Tuple
from Bio import Align, SeqIO
from Bio.Align import substitution_matrices
from Bio.Align.Applications import ClustalwCommandline, MuscleCommandline
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
import tempfile
import os


class AlignmentManager:
    """
    A class for managing sequence alignment operations.
    
    Supports multiple alignment algorithms including pairwise and multiple
    sequence alignment using Biopython's alignment tools.
    """
    
    def __init__(self, logger: Optional[logging.Logger] = None):
        """
        Initialize the AlignmentManager.
        
        Args:
            logger: Optional logger instance. If None, creates a default logger.
        """
        self.logger = logger or logging.getLogger(__name__)
        self.aligner = Align.PairwiseAligner()
        self.substitution_matrices = substitution_matrices
        
    def set_aligner_params(self, match_score: float = 2.0, 
                          mismatch_score: float = -1.0,
                          gap_open_score: float = -0.5,
                          gap_extend_score: float = -0.1,
                          substitution_matrix: str = 'BLOSUM62') -> None:
        """
        Set parameters for the pairwise aligner.
        
        Args:
            match_score: Score for matching characters
            mismatch_score: Score for mismatching characters
            gap_open_score: Score for opening a gap
            gap_extend_score: Score for extending a gap
            substitution_matrix: Name of substitution matrix to use
        """
        self.aligner.match_score = match_score
        self.aligner.mismatch_score = mismatch_score
        self.aligner.open_gap_score = gap_open_score
        self.aligner.extend_gap_score = gap_extend_score
        
        if substitution_matrix in self.substitution_matrices:
            self.aligner.substitution_matrix = self.substitution_matrices[substitution_matrix]
            self.logger.info(f"Using substitution matrix: {substitution_matrix}")
        else:
            self.logger.warning(f"Substitution matrix {substitution_matrix} not found, "
                              f"using default parameters")
    
    def pairwise_align(self, seq1: SeqRecord, seq2: SeqRecord) -> Tuple[float, str, str]:
        """
        Perform pairwise alignment between two sequences.
        
        Args:
            seq1: First sequence record
            seq2: Second sequence record
            
        Returns:
            Tuple containing (score, aligned_seq1, aligned_seq2)
        """
        try:
            alignments = self.aligner.align(seq1.seq, seq2.seq)
            best_alignment = alignments[0]
            
            aligned_seq1 = str(best_alignment.query)
            aligned_seq2 = str(best_alignment.target)
            score = best_alignment.score
            
            self.logger.debug(f"Pairwise alignment score: {score}")
            return score, aligned_seq1, aligned_seq2
            
        except Exception as e:
            self.logger.error(f"Error in pairwise alignment: {str(e)}")
            raise
    
    def multiple_sequence_align(self, sequences: List[SeqRecord], 
                               method: str = 'muscle') -> MultipleSeqAlignment:
        """
        Perform multiple sequence alignment.
        
        Args:
            sequences: List of SeqRecord objects to align
            method: Alignment method ('muscle', 'clustalw', or 'simple')
            
        Returns:
            MultipleSeqAlignment object
            
        Raises:
            ValueError: If method is not supported
            Exception: For alignment errors
        """
        if len(sequences) < 2:
            raise ValueError("At least 2 sequences are required for alignment")
        
        self.logger.info(f"Performing multiple sequence alignment using {method} "
                        f"on {len(sequences)} sequences")
        
        if method.lower() == 'muscle':
            return self._muscle_alignment(sequences)
        elif method.lower() == 'clustalw':
            return self._clustalw_alignment(sequences)
        elif method.lower() == 'simple':
            return self._simple_alignment(sequences)
        else:
            raise ValueError(f"Unsupported alignment method: {method}")
    
    def _muscle_alignment(self, sequences: List[SeqRecord]) -> MultipleSeqAlignment:
        """
        Perform alignment using MUSCLE.
        
        Args:
            sequences: List of SeqRecord objects
            
        Returns:
            MultipleSeqAlignment object
        """
        try:
            # Create temporary files for input and output
            with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as input_file:
                SeqIO.write(sequences, input_file.name, 'fasta')
                input_path = input_file.name
            
            with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as output_file:
                output_path = output_file.name
            
            # Run MUSCLE
            muscle_cline = MuscleCommandline(input=input_path, out=output_path)
            stdout, stderr = muscle_cline()
            
            # Read the alignment result
            alignment = Align.read(output_path, 'fasta')
            
            # Clean up temporary files
            os.unlink(input_path)
            os.unlink(output_path)
            
            self.logger.info("MUSCLE alignment completed successfully")
            return alignment
            
        except Exception as e:
            self.logger.error(f"MUSCLE alignment failed: {str(e)}")
            raise
    
    def _clustalw_alignment(self, sequences: List[SeqRecord]) -> MultipleSeqAlignment:
        """
        Perform alignment using ClustalW.
        
        Args:
            sequences: List of SeqRecord objects
            
        Returns:
            MultipleSeqAlignment object
        """
        try:
            # Create temporary files for input and output
            with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as input_file:
                SeqIO.write(sequences, input_file.name, 'fasta')
                input_path = input_file.name
            
            with tempfile.NamedTemporaryFile(mode='w', suffix='.aln', delete=False) as output_file:
                output_path = output_file.name
            
            # Run ClustalW
            clustalw_cline = ClustalwCommandline("clustalw2", infile=input_path, outfile=output_path)
            stdout, stderr = clustalw_cline()
            
            # Read the alignment result
            alignment = Align.read(output_path, 'clustal')
            
            # Clean up temporary files
            os.unlink(input_path)
            os.unlink(output_path)
            
            self.logger.info("ClustalW alignment completed successfully")
            return alignment
            
        except Exception as e:
            self.logger.error(f"ClustalW alignment failed: {str(e)}")
            raise
    
    def _simple_alignment(self, sequences: List[SeqRecord]) -> MultipleSeqAlignment:
        """
        Perform a simple progressive alignment.
        
        Args:
            sequences: List of SeqRecord objects
            
        Returns:
            MultipleSeqAlignment object
        """
        try:
            # Simple progressive alignment
            if len(sequences) == 2:
                score, aligned_seq1, aligned_seq2 = self.pairwise_align(sequences[0], sequences[1])
                
                # Create aligned sequences with gaps
                aligned_records = [
                    SeqRecord(Seq(aligned_seq1), id=sequences[0].id, description=sequences[0].description),
                    SeqRecord(Seq(aligned_seq2), id=sequences[1].id, description=sequences[1].description)
                ]
                
                return MultipleSeqAlignment(aligned_records)
            else:
                # For more than 2 sequences, fall back to a basic alignment
                # This creates a simple alignment by padding shorter sequences
                self.logger.warning("Simple alignment for >2 sequences using basic padding method.")
                
                max_length = max(len(seq.seq) for seq in sequences)
                aligned_records = []
                
                for seq_record in sequences:
                    # Pad sequence to max length
                    padded_seq = str(seq_record.seq).ljust(max_length, '-')
                    aligned_record = SeqRecord(
                        Seq(padded_seq), 
                        id=seq_record.id, 
                        description=seq_record.description
                    )
                    aligned_records.append(aligned_record)
                
                return MultipleSeqAlignment(aligned_records)
                
        except Exception as e:
            self.logger.error(f"Simple alignment failed: {str(e)}")
            raise
    
    def calculate_alignment_statistics(self, alignment: MultipleSeqAlignment) -> Dict[str, any]:
        """
        Calculate statistics for a multiple sequence alignment.
        
        Args:
            alignment: MultipleSeqAlignment object
            
        Returns:
            Dictionary containing alignment statistics
        """
        if not alignment:
            return {}
        
        alignment_length = alignment.get_alignment_length()
        num_sequences = len(alignment)
        
        # Calculate conservation
        conserved_positions = 0
        gap_positions = 0
        
        for i in range(alignment_length):
            column = alignment[:, i]
            unique_chars = set(char.upper() for char in column if char != '-')
            
            if len(unique_chars) <= 1:
                conserved_positions += 1
            elif '-' in column:
                gap_positions += 1
        
        conservation_percentage = (conserved_positions / alignment_length) * 100
        gap_percentage = (gap_positions / alignment_length) * 100
        
        stats = {
            "num_sequences": num_sequences,
            "alignment_length": alignment_length,
            "conserved_positions": conserved_positions,
            "conservation_percentage": conservation_percentage,
            "gap_positions": gap_positions,
            "gap_percentage": gap_percentage
        }
        
        self.logger.info(f"Alignment statistics: {num_sequences} sequences, "
                        f"length: {alignment_length}, "
                        f"conservation: {conservation_percentage:.1f}%")
        
        return stats
    
    def save_alignment(self, alignment: MultipleSeqAlignment, 
                      output_path: str, format: str = 'fasta') -> None:
        """
        Save alignment to a file.
        
        Args:
            alignment: MultipleSeqAlignment object to save
            output_path: Path where to save the alignment
            format: Output format ('fasta', 'clustal', 'phylip', etc.)
        """
        try:
            # Use SeqIO.write for MultipleSeqAlignment objects
            SeqIO.write(alignment, output_path, format)
            self.logger.info(f"Alignment saved to {output_path} in {format} format")
        except Exception as e:
            self.logger.error(f"Error saving alignment: {str(e)}")
            raise
