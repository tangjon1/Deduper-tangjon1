# Deduper Pseudocode

(Deduper Part 1)

## Define the problem

PCR duplicates from amplification bias can cause overrepresentation of specific sequences from a sequencing run. Using the unique molecular indexes of the adapter regions of the reads, PCR duplicates can be bioinformatically removed from the data to account for this bias. Furthermore, by using aligned reads, additional information from the alignment file—such as chromosome number, strandedness, and 5' start position of the read—can also be utilized to identify duplicates, drastically reducing computation time.

## Psuedocode

```
Load the valid UMIs into a set.
Instatiate a set to track the already processed unique reads.

Open the output file for writing.
Open the input SAM file for reading, store header reference lengths in a dictionary and write them to output (if available).
    For each aligned read record:
        Split the record on the tab delimiter.
        Get the RNAME of the read (str).
        Get the strandedness of the read (bool).
        [Get the UMI from the QNAME of the read (str).] <- Defined in a function
        [Get the 5' start position of the read from its POS, CIGAR string, and reference length (int).] <- Defined in a function
        [Get the true length of the read from the CIGAR string (int).] <- Defined in a function
        Construct the above elements into a tuple.
        Check if the tuple is in the set.
        If the tuple is not in the set:
            Write the record to the output file.
            Add the tuple to the set.

Close the output file.
```

## High Level Functions

### `get_umi()`

#### Description

```python
def get_umi(qname: str, valid_umis: set[str]) -> str | None:
    '''Return the UMI from the QNAME.
    Attempt to apply error correcting by Hamming distance = 1 to invalid UMIs.
    Returns None if the UMI fails error correction.'''
    return umi
```

#### Example

```python
>>> get_umi("NS500451:154:HWKTMBGXX:1:11101:15364:1139:GAACAGGT", valid_umis)
GAACAGGT
>>> get_umi("NS500451:154:HWKTMBGXX:1:11101:15364:1139:GAACAGG", valid_umis)
None
```

### `get_read_len()`

#### Description

```python
def get_read_len(cigar: str, true_len: bool=False) -> int:
    '''Parse the length of the read from a CIGAR string. Get the length of aligned reads with true_len = False, and the actual length with true_len = True.
    (True length is used to distinguish duplicates)'''
    return read_len
```

#### Example

```python
>>> get_read_len("10M")
10
>>> get_read_len("5M1000I5M")
10
>>> get_read_len("5M1000I5M", True)
1010
```

### `get_five_prime_start_pos()`

#### Description

```python
def get_five_prime_start_pos(record_pos: int, cigar: str, ref_len: int, strand: bool) -> int:
    '''Calculate the 5' start position of a stranded read with regard to the reference.
    Derive read_len from cigar -> (ref_len - (read_len + record_pos)) [for negative strands]'''
    return five_prime_start_pos
```

#### Example

```python
>>> get_five_prime_start_pos(5, "10M", 20, True)
5
>>> get_five_prime_start_pos(5, "10M", 20, False)
15
```
