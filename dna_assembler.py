#!/usr/bin/env python

"""
Procedures to assemble and analyze a small DNA sequence. Sequence reads are inputted
via a FASTQ file format, assembled into a single DNA sequence, then outputted as a
FASTA format file. The assembled reads are then submitted to the NCBI Blast web
service to determine the species that the sequence is from. Created for the Brandeis
University GPS course RSEG176 Cloud Computing, Spring 2016, Group B.

Aspects and portions of this code were derived from Ben Langmead's work in computational
genomics at Johns Hopkins University. This code and other information is available at
Ben's github repository at: https://github.com/BenLangmead/comp-genomics-class

Information on the NCBI BLAST program is available at:
http://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome
"""

# Utilizes the itertools library and the BioPython library
import itertools
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

__author__ = "Stephen S. Montanus"
__copyright__ = ""
__credits__ = ["Stephen S. Montanus"]
__license__ = "GPL"
__version__ = "0.5"
__maintainer__ = "Stephen S. Montanus"
__email__ = "smontanus@brandeis.edu"
__status__ = "Prototype"

ALLOWED_EXTENSIONS = set(['fq', 'FQ', 'fastq', 'FASTQ'])  # Check selected file for correct extension.
MAX_SEQLEN = 300 # Constraint for testing.
MAX_READS = 500 # Constraint for testing500.


def allowed_file(filename):
    """
    Checks the uploaded file name for correct extension and splits out the
    file name without the extension for naming FASTA files.

    Parameters:
        filename : Name of the FASTQ file, a string.

    Returns:
        A string representing the file name the FASTQ file without the extension.
    """
    return '.' in filename and \
           filename.rsplit('.', 1)[1] in ALLOWED_EXTENSIONS


def validfastq(fqfilename):
    """
    Validates the uploaded file as a FASTQ file, checks the data for correctness, and
    verifies that the file contents meets the constraints of the application for maximum
    sequence length and maximum number of reads.

    Parameters:
        fqfilename : Name of the FASTQ file, a string.

    Returns:
        A set of boolean values representing file contents, sequence length, and number of reads.
    """
    file_valid = True
    seqs_valid = True
    reads_valid = True
    num_reads = 0
    nucleotides = ["A", "C", "G", "T", "U", "N", "a", "c", "g", "t", "u", "n"]
    with open(fqfilename) as fqfilehandle:
            while True:
                fqfilehandle.readline()  # Skip the name line.
                seq = fqfilehandle.readline().rstrip()  # Read the base sequence.
                fqfilehandle.readline()  # Skip the placeholder line.
                fqfilehandle.readline().rstrip()  # Skip the base quality line.
                for n in seq:
                    if n not in nucleotides:
                        file_valid = False
                        return file_valid, seqs_valid, reads_valid
                if len(seq) > MAX_SEQLEN:
                    seqs_valid = False
                    return file_valid, seqs_valid, reads_valid
                num_reads += 1
                if num_reads > MAX_READS:
                    reads_valid = False
                    return file_valid, seqs_valid, reads_valid
                if len(seq) == 0:  # If the reads sequence is empty break from loop, EOF.
                    break
    return file_valid, seqs_valid, reads_valid


def readfastq(fqfilename):
    """
    Opens a FASTQ format file of DNA sequence reads and stores only the base
    sequence lines in a list for processing.

    Parameters:
        fqfilename : Name of the FASTQ file, a string.

    Returns:
        A list representing the individual reads from the FASTQ file as items.
    """
    reads = []  # Initialize the empty reads list.
    with open(fqfilename) as fqfilehandle:
        while True:
            fqfilehandle.readline()  # Skip the name line.
            seq = fqfilehandle.readline().rstrip()  # Read the base sequence.
            fqfilehandle.readline()  # Skip the placeholder line.
            fqfilehandle.readline().rstrip()  # Skip the base quality line.
            if len(seq) == 0:  # If the reads sequence is empty break from loop, EOF.
                break
            reads.append(seq)  # Append individual read sequence to the list as an item.
    return reads


def writefasta(fafilename, seq):
    """
    Writes a DNA sequence to a file in the FASTA format.

    Parameters:
        fafilename : Name of the FASTA file, a string.
        seq : DNA sequence to save, a string.

    Returns:
        Null.
    """
    fafilehandle = open(fafilename, 'w')  # Open file for writing in FASTA format.
    fafilehandle.write(">sequence_unknown")  # Write the FASTA header.
    fafilehandle.write("\n")  # Newline.
    fafilehandle.write(seq)  # Write the sequence data.
    fafilehandle.write("\n")  # Newline.
    fafilehandle.close()  # Close the file.


def readfasta(fafilename):
    """
    Opens a FASTA format file of a DNA sequence and stores only the
    sequence lines to a string for processing.

    Parameters:
        fafilename : Name of the FASTA file, a string.

    Returns:
        A string representing the DNA sequence from the FASTA file.
    """
    sequence = ''  # Initialize empty string to hold sequence data.
    with open(fafilename, 'r') as fafilehandle: # Open FASTA file.
        for line in fafilehandle:
            if not line[0] == '>':  # Ignore header line with genome information.
                sequence += line.rstrip()  # Concatenate lines into single sequence.
    return sequence


def analyzereads(readlist):
    """
    Performs a quick analysis of FASTQ file reads.

    Parameters:
        readlist : FASTQ file reads, a list.

    Returns:
        An integer representing the number of reads and a float representing the
        average read length.
    """
    numreads = len(readlist)
    readlengths = 0
    for i in readlist:
        readlengths += len(i)
    avglength = float(readlengths/numreads)
    return numreads, avglength


def overlap(seq_a, seq_b, min_length=3):
    """
    Compares two DNA sequence reads to determine if the suffix of one overlaps the
    prefix of the other.

    Parameters:
        seq_a : A DNA sequence read, a string.
        seq_b : A DNA sequence read, a string.
        min_length : Minimum length of overlap to search for, an integer.

    Returns:
        An integer representing the length of suffix prefix overlap or 0 if none exists.
    """
    start = 0  # Start search at the left.
    while True:
        start = seq_a.find(seq_b[:min_length], start)  # Look for seq_b's suffix in seq_a.
        if start == -1:  # No more character occurrences to right, so no match.
            return 0
        # Found occurrence, so check for full suffix/prefix match.
        if seq_b.startswith(seq_a[start:]):
            return len(seq_a)-start
        start += 1  # Move just past previous match.


def pick_maximal_overlap(reads, k):
    """
    Return a pair of reads from the list with a maximal suffix/prefix overlap >= k.
    Returns overlap length 0 if there are no such overlaps.

    Parameters:
        reads : A list DNA sequence reads, a list of strings.
        k : Minimum length of overlap to search for, an integer.

    Returns:
        Two strings representing the reads with longest overlap length >= k and the
        length of the overlap or 0 if none exists.
    """
    reada, readb = None, None  # Initialize empty strings for the reads.
    best_olen = 0  # Initialize variable for longest overlap length.
    for seq_a, seq_b in itertools.permutations(reads, 2):
        olen = overlap(seq_a, seq_b, min_length=k)
        if olen > best_olen:
            reada, readb = seq_a, seq_b
            best_olen = olen
    return reada, readb, best_olen


def greedy_scs(reads, k):
    """
    Greedy shortest-common-superstring merge.
    Repeat until no edges (overlaps of length >= k)
    remain.

    Parameters:
        reads : A list of DNA sequence reads, a list of strings.
        k : Minimum length of overlap to search for, an integer.

    Returns:
        A string representing the assembled reads.
    """
    read_a, read_b, olen = pick_maximal_overlap(reads, k)
    i = 0
    while olen > 0:
        i += 1
        reads.remove(read_a)
        reads.remove(read_b)
        reads.append(read_a + read_b[olen:])
        read_a, read_b, olen = pick_maximal_overlap(reads, k)
    return ''.join(reads)


def searchBlast(fastafile, e):
    """
    Submits a FASTA DNA sequence file to the NCBI Blast web service
    to determine the genome that the sequence belongs to.

    Parameters:
        fastafile : A DNA sequence in the FASTA file format, a string.
        e : E value threshold, e values of the alignment must be below this value, a float.

    Returns:
        A list with each item a list of alignment elements.
    """
    results = []
    fasta_string = open(fastafile).read()
    result_handle = NCBIWWW.qblast("blastn","nt", fasta_string)
    blast_record = NCBIXML.read(result_handle)
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            alignment_record = []
            if hsp.expect < e:
                alignment_record.append(alignment.title)
                alignment_record.append(alignment.length)
                alignment_record.append(hsp.expect)
                alignment_record.append(hsp.query)
                alignment_record.append(hsp.match)
                alignment_record.append(hsp.sbjct)
                results.append(alignment_record)
    return results