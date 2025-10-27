#!/usr/bin/env python

import argparse
import gzip
from typing import Any, IO, Callable
from io import TextIOWrapper

# For median function
from statistics import median

# Establish the CIGAR operations that consume reference
REF_CONS_CIG_OPS: set = {'M', 'D', 'N'}


def open_read(file_path: str) -> IO[Any] | TextIOWrapper:
    '''Takes an input file path and returns the appropriate open() function based on whether or not the file is compressed (ends with '.gz') or not.'''
    if file_path.endswith(".gz"):
        return gzip.open(file_path, 'rt')
    return open(file_path, 'r')


def get_args() -> argparse.Namespace:
    parser: argparse.ArgumentParser = argparse.ArgumentParser(description="Reference based deduplication to remove PCR duplicates from SAM files introduced by amplification bias")
    
    # Standard deduplication arguments 
    parser.add_argument('-f', '--file', help="The input SAM file, with uniquely mapped reads to deduplicate, sorted by reference", required=True)
    parser.add_argument('-o', '--outfile', help="Output path for resulting deduplicated SAM file", required=True)
    parser.add_argument('-u', '--umi', help="File containing UMI sequences, separated by newlines, to validate against")

    # Unimplemented arguments
    parser.add_argument('-c', '--correction', help="[UNIMPLEMENTED] The Hamming distance to use as an error correction threshold for UMI validation, where higher values are more lenient and lower values are more strict (with 0 requiring an exact match)", type=int, default=0)
    parser.add_argument('-r', '--represent', help="Determine how a read is selected for representation among duplicates", choices=['first', 'last', 'quality_mean', 'quality_median'], default='first')
    
    parser.add_argument('-s', '--force-single', help="Force the program to interpret reads as single-end reads", action="store_true")
    return parser.parse_args()


def validate_args(args: argparse.Namespace) -> None:
    '''TODO'''
    if args.correction != 0:
        raise argparse.ArgumentError(args.correction, "Error correction is unsupported; check the 'extra-features' branch!")
    if args.represent != 'first':
        raise argparse.ArgumentError(args.represent, "Selecting representatives other than the first encountered read is unsupported; check the 'extra-features' branch!")
    if args.paired:
        raise argparse.ArgumentError(args.paired, "Paired-end read deduplication is unsupported; check the 'extra-features' branch!")
    if args.umi is None:
        raise argparse.ArgumentError(args.umi, "Randomers are not supported and UMIs must be provided; check the 'extra-features' branch!")

    if args.umi is None and args.correction != 0:
            raise argparse.ArgumentError(args.correction, "Error correction can only be applied using known UMIs")
    if args.correction and args.correction < 0:
        raise argparse.ArgumentError(args.correction, "Error correction threshold cannot be less than 0")


def qscore_none(*_) -> None:
    '''TODO'''
    return None


def qscore_mean(qscore_seq: str, phred_offset: int=33) -> float | None:
    '''TODO'''
    if qscore_seq == '*':
        return None
    sum: int = 0
    for char in qscore_seq:
        sum += ord(char) - phred_offset
    return sum / len(qscore_seq)


def qscore_median(qscore_seq: str, phred_offset: int=33) -> float | None:
    '''TODO'''
    if qscore_seq == '*':
        return None
    qscore_vals: list[int] = []
    for char in qscore_seq:
        qscore_vals.append(ord(char) - phred_offset)
    return median(qscore_vals)


def import_cigar_string(cigar_string: str) -> tuple[tuple[str, int], ...]:
    '''
    Convert a CIGAR string into a Python-parsable format:
    A tuple that contains tuples of pairs of values, where the first value is the operation and the second value is its associated number of bases
    Order is preserved from the original string. 
    '''
    cigar_seq: list[tuple[str, int]] = list()
    cigar_int: str = ""
    for char in cigar_string:
        if char.isdigit():
            cigar_int += char
        else:
            cigar_seq.append((char, int(cigar_int)))
            cigar_int = ""
    return tuple(cigar_seq)


def get_five_prime_softclip(cigar: tuple[tuple[str, int], ...], reverse_strand: bool) -> int:
    '''Get the number of 5' soft clipped bases from an imported CIGAR, given strandedness'''
    # Ignore hard clipping
    hardclip_start: bool = cigar[0][0] == 'H'
    hardclip_end: bool = cigar[-1][0] == 'H'
    if hardclip_start or hardclip_end:
        cigar = cigar[hardclip_start:len(cigar) - hardclip_end]
    
    # Get the soft clipping value from the appropriate side of the cigar
    if cigar[-reverse_strand][0] == 'S':
        return cigar[-reverse_strand][1]

    # No soft clipping detected, return 0
    return 0
    

def get_file_max_five_prime_adjustment(file: str) -> int:
    '''Get the maximum difference between POS and 5' POS that appears among the reads of a SAM file'''
    max_adj_pos_diff: int = 0
    with open_read(file) as in_file:
        for line in in_file:
            # Skip the header lines
            if line[0] == '@':
                continue

            sam_record_split: list[str] = line.strip('\n').split('\t')
            flag: int = int(sam_record_split[1])
            reverse_strand: bool = flag & 16 == 16
            pos = int(sam_record_split[3])

            # Parse the CIGAR strings
            cigar_string: str = line.strip('\n').split('\t')[5]
            if cigar_string == '*':
                continue
            else:
                cigar = import_cigar_string(cigar_string)
                max_adj_pos_diff = max(max_adj_pos_diff, abs(get_five_prime_start(pos, cigar, reverse_strand) - pos))
    return max_adj_pos_diff


def get_cigar_ref_cons_length(cigar: tuple[tuple[str, int], ...]) -> int:
    '''Get the sum of all reference-consuming operations of a cigar'''
    ref_cons_length: int = 0
    for operation in cigar:
        if operation[0] in REF_CONS_CIG_OPS:
            ref_cons_length += operation[1]
    return ref_cons_length


def build_valid_umi_set(umi_filepath: str) -> set[str]:
    '''Given a file where each line represents a valid UMI, construct and return a set containing those UMIs'''
    valid_umis = set()
    with open_read(umi_filepath) as umi_file:
        for line in umi_file:
            valid_umis.add(line.strip('\n'))
    return valid_umis


def get_nearest_valid_umi(umi: str, error_correction: int, valid_umis: set[str] | None) -> str | None:
    '''TODO'''
    raise NotImplementedError


def get_umi(qname: str, error_correction: int, valid_umis: set[str] | None, second_read: bool | None=None) -> str | None:
    '''
    Extract the UMI from a QNAME field and validate it against the set of valid UMIs, returning the UMI if it is valid.
    For single-end reads, set second_read to None. For paired-end reads, set second_read to False/True for the R1/R2 reads, respectively.
    If the UMI could not be validated, returns None.
    '''
    umi_index = -1
    if second_read is not None:
         umi_index -= (not second_read)
    
    umi: str = qname.split(':')[umi_index]

    if valid_umis:
        if umi in valid_umis:
            return umi
        elif error_correction > 0:
            return get_nearest_valid_umi(umi, error_correction, valid_umis)
        else:
            return None
    else:
        return umi


def get_five_prime_start(pos: int, cigar: tuple[tuple[str, int], ...], reverse_strand: bool) -> int:
    '''TODO'''
    five_prime_softclip: int = get_five_prime_softclip(cigar, reverse_strand)
    if not reverse_strand:
        return pos - five_prime_softclip
    else:
        rev_pos_offset: int = get_cigar_ref_cons_length(cigar)
        return pos + rev_pos_offset + five_prime_softclip - 1
    

def process_sam_record(sam_record: str, error_correction: int, valid_umis: set[str] | None, qscore_method: Callable[[str], int | float | None]) -> tuple[tuple[str, int, bool, str, bool | None], str, int | float | None] | None:
    '''
    Process the necessary fields of a SAM read record into a tuple.
    Returns None if the record is determined to be invalid for any reason.
    '''
    # Split the record on the tab delimiter, parse bitwise flag
    sam_record_split: list[str] = sam_record.split('\t')
    flag: int = int(sam_record_split[1])

    # Check the UMI in the read name (QNAME), if get_umi() returns None
    qname: str = sam_record_split[0]
    paired: bool = flag & 1 == 1
    if paired:
        second_read: bool | None = flag & 128 == 128
    else:
        second_read: bool | None = None

    umi = get_umi(qname, error_correction, valid_umis, second_read)
    if umi is None:
        return None
    
    # Get the reference name (RNAME) from the record
    rname: str = sam_record_split[2]

    # Determine strandedness from the bitflag (FLAG)
    reverse_strand: bool = flag & 16 == 16

    # Determine the 5' start position from left-based position (POS), the CIGAR string (CIGAR), and strandedness
    pos: int = int(sam_record_split[3])
    cigar: tuple[tuple[str, int], ...] = import_cigar_string(sam_record_split[5])
    five_prime_start: int = get_five_prime_start(pos, cigar, reverse_strand)

    # Calculate quality score metric with the provided method
    qscore: int | float | None = qscore_method(sam_record_split[10])

    # QNAME returned for paired read mate identification purposes
    return (rname, five_prime_start, reverse_strand, umi, second_read), qname, qscore


def main(args: argparse.Namespace) -> None:
    '''TODO'''
    # Get the maximum leftmost soft clipping value
    max_adj_pos_diff: int = get_file_max_five_prime_adjustment(args.file)
    print(f"max_adj_pos_diff\t{max_adj_pos_diff}")

    # Build a set of valid UMIs
    valid_umis: set[str] = build_valid_umi_set(args.umi)

    # Track the unique encountered reads for duplication in a dictionary;
    # Keys represent the attributes that define uniqueness as a tuple
    # Values represent the a tuple containing:
    #   - quality score metric (if applicable, otherwise None)
    #   - the sam record line
    unique_reads: dict[tuple[tuple[str, int, bool, str, bool | None], ...], tuple[int | float | None, str]] = dict()

    # For paired end reads, we need to track the mates that have already been parsed;
    # track the key QNAME (str) the R1/R2 status (bool), quality metric (numeric),
    # and entire record of the read (str)
    unique_mate_pair_qname: dict[str, dict[tuple[str, int, bool, str, bool | None], tuple[int | float | None, str]]] = dict()

    # Establish the quality score method from the user-provided arguments
    qscore_method = qscore_none
    if args.represent == 'quality_mean':
        qscore_method = qscore_mean
    elif args.represent == 'quality_median':
        qscore_method = qscore_median
        
    # Initialize the location (RNAME and adj. POS) to None
    # adj. POS is POS offset by the maximum soft clipping value, for a "sliding window" (better memory efficiency)
    current_reference = None
    current_adj_pos = None

    # Track various stats about the deduplicated file
    stats: dict[str, int] = {
        'headers': 0,
        'unique_reads': 0,
        'invalid_umi_records': 0,
        'duplicates_removed': 0
    }

    # Track the number of reads per RNAME
    rname_count: dict[str, int] = dict()

    with (
        open_read(args.file) as in_file,
        open(args.outfile, 'w') as out_file
    ):
        for line in in_file:
            # Write header lines to output file
            if line[0] == '@':
                stats['headers'] += 1
                out_file.write(line)
                continue

            # Process each SAM record for read data
            sam_record: str = line.strip('\n')
            processed_sam_record = process_sam_record(sam_record, args.correction, valid_umis, qscore_method)
            if processed_sam_record is None:
                # If processing returns None, the read was invalid -- continue to next read
                stats['invalid_umi_records'] += 1
                continue
            
            # Expand the record tuple into the identifying read data, QNAME, and quality score metric
            sam_read, qname, qscore = processed_sam_record

            # Main duplication checking logic
            sam_read_value = qscore, sam_record

            read_to_check: tuple[tuple[str, int, bool, str, bool | None], ...] = (sam_read,)

            # For paired reads, check if the mate has already been found
            # If the mate has been found, construct the new read_to_check from the pair
            paired = sam_read[4] is not None
            if paired and not args.force_single:
                raise NotImplementedError("Paired end reads not yet supported")
                if qname in unique_mate_pair_qname:      
                    paired_sam_read: list[tuple[str, int, bool, str, bool | None] | None] = [None, None]
                    # This loop should only iterate once
                    if len(unique_mate_pair_qname[qname]) != 1:
                        raise ValueError("This dictionary should contain only 1 element")
                    for sam_read_key in unique_mate_pair_qname[qname]:
                        if sam_read_key[4] is None:
                            raise ValueError("Paired read does not have appropriate bool set")

                        paired_sam_read[sam_read_key[4]] = sam_read_key
                        paired_sam_read[not sam_read_key[4]] = sam_read

                        # This is mainly to appease type annotations
                        if paired_sam_read[0] is None or paired_sam_read[1] is None:
                            raise ValueError("Read tuples should not be None")
                        
                        # Construct the new read_to_check
                        read_to_check = paired_sam_read[0], paired_sam_read[1]
                
                        # Add the qscores and concatenate the read record strings
                        mate_qscore = unique_mate_pair_qname[qname][sam_read_key][0]
                        new_qscore: int | float | None = None
                        if qscore is not None and mate_qscore is not None:
                            new_qscore = qscore + mate_qscore
                        
                        if sam_read_key[4]:
                            new_sam_record: str = sam_record + '\n' + unique_mate_pair_qname[qname][sam_read_key][1]
                        else:
                            new_sam_record: str = unique_mate_pair_qname[qname][sam_read_key][1] + '\n' + sam_record

                        # Get the new sam read values from those of those of each read in the pair
                        sam_read_value = new_qscore, new_sam_record

                    # Free up memory since we do not expect that QNAME elsewhere in file after finding the mate
                    del unique_mate_pair_qname[qname]
                else:
                    unique_mate_pair_qname[qname] = {sam_read: (qscore, sam_record)}

            # Check that the read is unique by searching for it in the unique_reads dictionary
            if read_to_check in unique_reads:
                # Track the number of discarded duplicates
                stats['duplicates_removed'] += 1 + paired

                # Check if the existing read needs to be updated according to the representation type...
                # ...for 'last' represent
                if args.represent == 'last':
                    unique_reads[read_to_check] = sam_read_value

                # ...for highest quality represent
                if 'quality' in args.represent:
                    prev_qscore = unique_reads[read_to_check][0]
                    # Error if either qscore is None
                    if prev_qscore is None or qscore is None:
                        raise ValueError("Quality score metric should not be None for this representation type")
                    if qscore > prev_qscore:
                        unique_reads[read_to_check] = sam_read_value
           
            else:
                # If the read is not already found, add it to the dictionary
                unique_reads[read_to_check] = sam_read_value
            
            # Update the location
            current_reference = sam_read[0]
            current_window_min: int = sam_read[1] - (max_adj_pos_diff * 2)

            # Write the out-of-window SAM records to the file
            for key in list(unique_reads.keys()):
                key_reference = key[0][0]
                key_adj_pos = key[0][1]

                if current_reference != key_reference or key_adj_pos < current_window_min:
                    out_file.write(unique_reads[key][1] + '\n')
                    del unique_reads[key]

                    # Track the number of total reads written, and per chromosome
                    stats['unique_reads'] += len(key)
                    rname_count.setdefault(key_reference, 0)
                    rname_count[key_reference] += len(key)

                # Since dictionaries are ordered, assume that when we encounter a read within the window,
                # all remaining reads should also be retained in the dictionary
                elif current_reference == key_reference and key_adj_pos >= current_window_min:
                    break
    
        # Dump remaining unique reads to output file
        for key in list(unique_reads.keys()):
            key_reference = key[0][0]
            key_adj_pos = key[0][1]
            out_file.write(unique_reads[key][1] + '\n')
            stats['unique_reads'] += len(key)
            rname_count.setdefault(key_reference, 0)
            rname_count[key_reference] += len(key)

    
    # Print stats to standard out
    for stat in stats:
        print(f"{stat}\t{stats[stat]}")

    for rname in rname_count:
        print(f"{rname}\t{rname_count[rname]}")

            
if __name__ == '__main__':
    args = get_args()
    main(args)
