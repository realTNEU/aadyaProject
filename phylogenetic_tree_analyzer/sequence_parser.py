"""
Sequence Parser Module

This module handles parsing of biological sequence data from various formats
including FASTA and GenBank files.
"""

import os
import logging
from typing import List, Dict, Optional, Union
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO.FastaIO import SimpleFastaParser


class SequenceParser:
    """
    A class for parsing biological sequence data from various file formats.
    
    Supports FASTA and GenBank formats with error handling and validation.
    """
    
    def __init__(self, logger: Optional[logging.Logger] = None):
        """
        Initialize the SequenceParser.
        
        Args:
            logger: Optional logger instance. If None, creates a default logger.
        """
        self.logger = logger or logging.getLogger(__name__)
        self.supported_formats = ['fasta', 'fa', 'genbank', 'gb', 'gbk']
        
    def parse_file(self, file_path: str) -> List[SeqRecord]:
        """
        Parse a sequence file and return a list of SeqRecord objects.
        
        Args:
            file_path: Path to the sequence file
            
        Returns:
            List of SeqRecord objects
            
        Raises:
            FileNotFoundError: If the file doesn't exist
            ValueError: If the file format is not supported
            Exception: For other parsing errors
        """
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"File not found: {file_path}")
            
        file_extension = os.path.splitext(file_path)[1].lower().lstrip('.')
        
        if file_extension not in self.supported_formats:
            raise ValueError(f"Unsupported file format: {file_extension}. "
                           f"Supported formats: {', '.join(self.supported_formats)}")
        
        try:
            self.logger.info(f"Parsing sequence file: {file_path}")
            sequences = list(SeqIO.parse(file_path, file_extension))
            
            if not sequences:
                self.logger.warning(f"No sequences found in file: {file_path}")
                return []
                
            self.logger.info(f"Successfully parsed {len(sequences)} sequences from {file_path}")
            return sequences
            
        except Exception as e:
            self.logger.error(f"Error parsing file {file_path}: {str(e)}")
            raise
    
    def parse_fasta(self, file_path: str) -> List[SeqRecord]:
        """
        Parse a FASTA file specifically.
        
        Args:
            file_path: Path to the FASTA file
            
        Returns:
            List of SeqRecord objects
        """
        return self.parse_file(file_path)
    
    def parse_genbank(self, file_path: str) -> List[SeqRecord]:
        """
        Parse a GenBank file specifically.
        
        Args:
            file_path: Path to the GenBank file
            
        Returns:
            List of SeqRecord objects
        """
        return self.parse_file(file_path)
    
    def validate_sequences(self, sequences: List[SeqRecord], 
                          min_length: int = 10) -> List[SeqRecord]:
        """
        Validate sequences and filter out invalid ones.
        
        Args:
            sequences: List of SeqRecord objects to validate
            min_length: Minimum sequence length required
            
        Returns:
            List of validated SeqRecord objects
        """
        validated_sequences = []
        
        for seq_record in sequences:
            if len(seq_record.seq) < min_length:
                self.logger.warning(f"Sequence {seq_record.id} is too short "
                                  f"({len(seq_record.seq)} < {min_length}), skipping")
                continue
                
            if not seq_record.seq:
                self.logger.warning(f"Sequence {seq_record.id} is empty, skipping")
                continue
                
            validated_sequences.append(seq_record)
        
        self.logger.info(f"Validated {len(validated_sequences)} out of {len(sequences)} sequences")
        return validated_sequences
    
    def get_sequence_info(self, sequences: List[SeqRecord]) -> Dict[str, any]:
        """
        Get summary information about a list of sequences.
        
        Args:
            sequences: List of SeqRecord objects
            
        Returns:
            Dictionary containing sequence statistics
        """
        if not sequences:
            return {"count": 0, "total_length": 0, "avg_length": 0}
        
        lengths = [len(seq.seq) for seq in sequences]
        total_length = sum(lengths)
        avg_length = total_length / len(lengths)
        
        info = {
            "count": len(sequences),
            "total_length": total_length,
            "avg_length": avg_length,
            "min_length": min(lengths),
            "max_length": max(lengths),
            "sequence_ids": [seq.id for seq in sequences]
        }
        
        self.logger.info(f"Sequence info: {info['count']} sequences, "
                        f"avg length: {avg_length:.1f}")
        
        return info
    
    def save_sequences(self, sequences: List[SeqRecord], 
                      output_path: str, format: str = 'fasta') -> None:
        """
        Save sequences to a file.
        
        Args:
            sequences: List of SeqRecord objects to save
            output_path: Path where to save the sequences
            format: Output format ('fasta', 'genbank', etc.)
        """
        try:
            SeqIO.write(sequences, output_path, format)
            self.logger.info(f"Saved {len(sequences)} sequences to {output_path}")
        except Exception as e:
            self.logger.error(f"Error saving sequences to {output_path}: {str(e)}")
            raise
