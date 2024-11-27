# Import required libraries for the operation
import os
import subprocess
import sys
import argparse
from collections import defaultdict
from pprint import pprint
import time
from datetime import datetime
import psutil

start_time = time.time()  # gets the current time

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('ref_genome', metavar='R', type=str, help='The path to the reference genome file')
args = parser.parse_args()

def log_ram_usage():
    # Get the process id
    pid = os.getpid()
    process = psutil.Process(pid)
    # Get the memory info
    mem_info = process.memory_info()
    rss_memory = mem_info.rss // 1024  # Convert to kB
    vsz_memory = mem_info.vms // 1024  # Convert to kB

    # Get the current timestamp
    timestamp = datetime.now()

    # Open the log file in append mode
    with open('ram_usage_log.txt', 'a') as log_fh:
        # Write the memory usage and timestamp to the log file
        log_fh.write(f"Timestamp: {timestamp}\n")
        log_fh.write(f"RSS Memory Usage: {rss_memory} kB\n")
        log_fh.write(f"VSZ Memory Usage: {vsz_memory} kB\n")
        log_fh.write("-----------------------\n")

def reconstruct_alignment(cigar_total_length, cigar_arrayref, readseq, readqual, refseq, refpos_start):
    """
    This function reconstructs alignment from a cigar string which denotes the alignment between a read sequence and reference sequence.

    Args:
        cigar_total_length: total length of the read
        cigar_arrayref: CIGAR string
        readseq: actual sequence of the read
        readqual: quality scores for the bases in the read
        refseq: reference genome sequence
        refpos_start: starting position of the read's alignment in the reference genome
    """
    try:
        output = []  # List to store the reconstructed alignment information

        # Convert sequences to lists with a None element at index 0 for 1-based indexing
        readseq = [None] + list(readseq)  
        readqual = [None] + list(readqual)  
        refseq = [None] + list(refseq)  

        cigar_op = [None]  # List to store CIGAR operations
        cigar_pos_offset = 1  # Offset for cigar operations

        cigar_readseq = [None]  # List to store read sequences based on CIGAR operations
        cigar_readqual = [None]  # List to store quality scores based on CIGAR operations
        cigar_refpos = [None]  # List to store reference positions based on CIGAR operations
        cigar_refseq = [None]  # List to store reference sequences based on CIGAR operations

        readbase_pos_offset = 0  # Offset for read base positions
        refpos_pos_offset = refpos_start  # Offset for reference positions

        # Iterate over each cigar operation to populate the alignment lists
        for cigarop in cigar_arrayref[0]:
            length = cigarop[0]  # Length of the cigar operation
            op = cigarop[1]  # Type of cigar operation

            # Populate alignment lists based on the cigar operation type and length
            for i in range(cigar_pos_offset, cigar_pos_offset + length):
                cigar_op.append(op)
                cigar_readseq.append('*' if op == "D" else readseq[i + readbase_pos_offset])
                cigar_readqual.append(' ' if op == "D" else readqual[i + readbase_pos_offset])
                cigar_refpos.append('*' if op == "I" or op == "S" else refpos_pos_offset + i - cigar_pos_offset)
                cigar_refseq.append('*' if op == "I" or op == "S" else refseq[refpos_pos_offset + i - cigar_pos_offset - refpos_start + 1])

            # Update offsets based on the current cigar operation
            if op == "D":
                readbase_pos_offset -= length
            if op != "I" and op != "S":
                refpos_pos_offset += length
            cigar_pos_offset += length

        # If debugging flag is set, print the alignment-related information
        if debug:
            print("\t".join(map(lambda x: x if x is not None else "_", cigar_refpos[1:cigar_total_length])))
            print("\t".join(map(lambda x: x if x is not None else "_", cigar_refseq[1:cigar_total_length])))
            print("\t".join(map(lambda x: x if x is not None else "_", cigar_op[1:cigar_total_length])))
            print("\t".join([cigar_op[i] if cigar_op[i] == "M" else (cigar_readseq[i] if cigar_readseq[i] == cigar_refseq[i] else "X") for i in range(1, cigar_total_length)]))
            print("\t".join(map(lambda x: x if x is not None else "_", cigar_readseq[1:cigar_total_length])))
            print("\t".join(map(lambda x: x if x is not None else "_", cigar_readqual[1:cigar_total_length])))
            print("\t".join([str(ord(x) - phredoffset) if x is not None else "-1" for x in cigar_readqual[1:cigar_total_length]]))

        start = -1  # Start position of an alignment window
        end = -1 * max_distance  # End position of an alignment window

        min_window_size = 2  # Minimum size requirement for an alignment window

        # Iterate over the cigar operations to identify alignment windows
        for i in range(1, cigar_total_length+1):
            if not ((cigar_op[i] == "M" and cigar_readseq[i] != cigar_refseq[i]) or cigar_op[i] == "I" or cigar_op[i] == "D"):
                continue

            # Check if the current position exceeds the maximum distance from the previous window
            if (i - end) > max_distance:
                # Check if a valid alignment window exists and meets the minimum window size requirement
                if start >= 0 and (end - start + 1) >= min_window_size:
                    filtered_seq = list(filter(lambda x: x != "*", cigar_refpos[start:end]))

                    if filtered_seq:  # Check if the filtered sequence is not empty
                        # Create an output entry with relevant alignment information
                        output.append([start, end, min(filtered_seq), max(filtered_seq),
                                       "".join(filter(lambda x: x != "*", cigar_refseq[start:end])),
                                       "".join(filter(lambda x: x != "*", cigar_readseq[start:end])),
                                       mean([ord(x) - phredoffset for x in filter(lambda x: x != " ", cigar_readqual[start:end])])])
                    else:
                        # Handle the case where the filtered sequence is empty
                        print(f"No valid positions found in region {start} to {end}.")

                start = i  # Update the start position for the new window
                end = i  # Update the end position for the new window
            else:
                if i > end:
                    end = i  # Update the end position of the current window

        # Check if the last alignment window meets the minimum window size requirement
        if start >= 0 and (end - start + 1) >= min_window_size:
            filtered_seq = list(filter(lambda x: x != "*", cigar_refpos[start:end]))

            if filtered_seq:  # Check if the filtered sequence is not empty
                # Create an output entry with relevant alignment information
                output.append([start, end, min(filtered_seq), max(filtered_seq),
                               "".join(filter(lambda x: x != "*", cigar_refseq[start:end])),
                               "".join(filter(lambda x: x != "*", cigar_readseq[start:end])),
                               mean([ord(x) - phredoffset for x in filter(lambda x: x != " ", cigar_readqual[start:end])])])
            else:
                # Handle the case where the filtered sequence is empty
                print(f"No valid positions found in region {start} to {end}.")

        return output  # Return the reconstructed alignment information as output

    except Exception as e:
        print("Error occurred in reconstruct_alignment function: ", str(e))
        sys.exit(1)  # Exit the script if an error occurs


def mean(array):
    # This function calculates the mean of the values in an array.
    # If the array is empty, it returns -1.
    if len(array) == 0:
        return -1
    else:
        return sum(array) / len(array)

def parse_cigar(cigar_string):
    """
    This function parses cigar strings that contain information about matching sequences, insertions, and deletions.
    It converts the cigar string into an array of [length, operation] pairs.
    The function uses regular expressions (re) to find all (length, operation) pairs in the cigar string.
    """
    try:
        import re
        cigar = []
        i = 0
        # Find all (length, operation) pairs in the cigar string using regular expression matching.
        matches = re.findall(r'([0-9]+)([MIDNSHPX=])', cigar_string)
        for match in matches:
            if match[1] == "N" or match[1] == "P":
                # Raise a ValueError if an unexpected CIGAR operation is encountered.
                raise ValueError("ERROR: Unexpected CIGAR operation " + match[1])
            if match[1] != "H":
                # Add the [length, operation] pair to the cigar array, except for "H" (hard clip) operations.
                cigar.append([int(match[0]), match[1]])
                i += int(match[0])
        return (i, cigar)
    except Exception as e:
        # Handle any exceptions that occur during the execution of the function.
        print("Error occurred in parse_cigar function: ", str(e))
        sys.exit(1)

# The following two functions (faidx and depth) cache results to save execution time in future calls.

faidxcache = {} # Creating a cache dictionary to store faidx results
faidxcache_count = 0 # Initializing a count to keep track of the number of cached results

def faidx(region):
    """
    This function uses samtools to fetch the sequence for a specific genomic region from a reference genome. It also implements caching to avoid repeated samtools calls for the same region.
    """
    try:
        global faidxcache, faidxcache_count # Giving the function access to the global cache variables
        # If the sequence for this region is already cached, return it
        if region in faidxcache:
            return faidxcache[region]
        newRegion = f"chr{region}" # Formatting region to work with samtools syntax
        output = subprocess.getoutput(f'{samtools} faidx {refseq} {newRegion}') # Calling samtools faidx command to fetch sequence
        outputlines = output.split("\n") # Splitting output into lines
        # Checking if the output is valid
        if outputlines[0].split(' ')[0] != ">" + newRegion:
            raise ValueError("Error: unexpected output from samtools faidx")
        seq = "".join(outputlines[1:]) # Joining the lines to form the sequence
        faidxcache[region] = seq # Storing sequence in cache
        faidxcache_count += 1 # Increasing cache count
        return seq # Returning sequence
    except Exception as e:
        print("Error occurred in faidx function: ", str(e)) # Printing error message and stopping execution if an error occurs
        sys.exit(1)

depthcache = {} # Creating a cache dictionary to store depth results
depthcache_count = 0 # Initializing a count to keep track of the number of cached results

def depth(region):
    """
    This function calculates the depth of coverage for a specific genomic region using samtools. It also implements caching to avoid repeated samtools calls for the same region.
    """
    try:
        global depthcache, depthcache_count # Giving the function access to the global cache variables
        # If the depth for this region is already cached, return it
        if region in depthcache:
            return depthcache[region]
        newRegion = f"chr{region}" # Formatting region to work with samtools syntax
        output = subprocess.getoutput(f'{samtools} depth -d {max_samtools_depth} -r {newRegion} -q {quality_threshold} {depth_bam}') # Calling samtools depth command to calculate depth
        outputlines = output.split("\n") # Splitting output into lines
        depth_values = [] # Creating list to store depth values
        for line in outputlines:
            fields = line.split("\t") # Splitting line into fields
            depth_values.append(int(fields[2])) # Appending depth value to list
        depthcache[region] = mean(depth_values) # Calculating mean depth and storing it in cache
        depthcache_count += 1 # Increasing cache count
        
        return depthcache[region] # Returning mean depth
    except Exception as e:
        print("Error occurred in depth function: ", str(e)) # Printing error message and stopping execution if an error occurs
        sys.exit(1)

# User is prompted for input to provide the paths to the BAM file and the reference genome.
try:
    bam_file = bam_file = sys.stdin # Reading the BAM file from standard input
    reference_genome = args.ref_genome # Getting the path to the reference genome from the command line arguments
except Exception as e:
    print("Error occurred while taking user inputs: ", str(e)) 
    sys.exit(1)

# Setting up some default parameters.
refseq = reference_genome  # Assigning the reference genome to refseq
samtools = "samtools"  # Setting the samtools command
rawoutput = 0  # Setting rawoutput to false. If true, outputs would not be prettified
debug = 0  # Setting debug mode to false. If true, additional debugging information would be printed
phredoffset = 33  # Setting the default PHRED quality score offset
quality_threshold = 20  # Setting the default quality threshold for considering base calls
min_depth = 50  # Setting the minimum depth of coverage required to consider a region
min_af = 0  # Setting the minimum allele frequency to consider a mutation
depth_bam = bam_file  # Assigning the bam file to depth_bam variable for depth calculation
skip_lowqual = 0  # If set to true, regions with quality score lower than quality_threshold will be skipped
skip_lowdepth = 0  # If set to true, regions with depth of coverage lower than min_depth will be skipped
skip_lowaf = 0  # If set to true, mutations with allele frequency lower than min_af will be skipped
max_samtools_depth = 500000  # overriding samtools depth default cap of 8000x depth to 500000x
max_distance = 5  # Setting the maximum distance between two mutations to consider them as linked
detected_mutations = {}  # Initializing a dictionary to store detected mutations

increment = 100000  # Setting the amount by which to increment the region size for each samtools call
line_number = 0  # Initializing a counter to keep track of the number of lines processed

ram_time = time.time() # Initialize the RAM check time

# Main loop that reads in lines of data from the input BAM file, processes each line, and outputs variant calls.
while True:
    try:
        line = bam_file.readline().strip() # Read a line from the BAM file
        line_number += 1 # Increment the line number
        # If the number of processed lines is a multiple of increment, print number of processed lines
        if line_number % increment == 0:
            print(f"Processed {line_number} lines", file=sys.stderr)

        # Check if a minute has passed and log RAM usage if so
        if time.time() - ram_time >= 60:
            log_ram_usage()
            # Reset the start time for the next minute
            ram_time = time.time()

        # If the line is empty or starts with '@' (header in BAM files), skip this iteration
        if len(line) == 0 or line[0] == "@":
            continue

        fields = line.split("\t") # Split the line into fields by tab character
        
        # If any of the key fields are missing ("*"), skip this iteration
        if fields[2] == "*" or fields[3] == "*" or fields[5] == "*":
            continue

        cigar_total_length, *cigar = parse_cigar(fields[5]) # Parse the CIGAR string to calculate the total length and extract individual operations
        rawrefseq = faidx(f"{fields[2]}:{fields[3]}-{int(fields[3]) + cigar_total_length - 1}") # Fetch the reference sequence corresponding to this read
        
        try:
            readseq = fields[9].upper() # Get the read sequence and convert it to uppercase
        except AttributeError:
            print(f"Missing sequence on line {line_number}. Skipping this line.") # If the sequence is missing, print a message and skip this line
            continue

        readqual = fields[10] # Get the quality scores of the read
        direction = "-" if int(fields[1]) & 0x0010 else "+" # Determine the direction of the read

        # If in debug mode, print some information about the read
        if debug:
            print("\t".join([fields[0], str(len(fields)), fields[5], readseq, readqual]))

        # Reconstruct the alignment to identify potential mutations
        candidate_mutations = reconstruct_alignment(cigar_total_length, cigar, readseq, readqual, rawrefseq.upper(), int(fields[3]))

        # Iterate over all candidate mutations
        for candidate_mutation_arrayref in candidate_mutations:
            # Get the lengths of the reference and variant alleles
            candidate_mutation_ref_len = len(candidate_mutation_arrayref[4])
            candidate_mutation_var_len = len(candidate_mutation_arrayref[5])
            
            if candidate_mutation_ref_len >= 1 and candidate_mutation_var_len == 0:
                # If it's a simple deletion, skip this iteration
                continue
            elif candidate_mutation_ref_len == 0 and candidate_mutation_var_len >= 1:
                # If it's a simple insertion, skip this iteration
                continue
            else:
                # If it's a complex indel, record it
                if rawoutput:
                    print("\t".join([fields[0]] + candidate_mutation_arrayref + [direction]))
                
                # Create a key to identify this mutation
                key = "|".join(str(item) for item in ([fields[2]] + candidate_mutation_arrayref[2:6]))
                
                # If this mutation is not yet in the detected_mutations dictionary, add it
                if key not in detected_mutations:
                    detected_mutations[key] = [0, 0, [], []]
                
                # Increment the count for the direction of this read
                idx = 1 if direction == "-" else 0
                detected_mutations[key][idx] += 1
                
                # Append the quality score of this mutation
                detected_mutations[key][2 + idx].append(candidate_mutation_arrayref[6])
    except Exception as e:
        # If an error occurred while processing the line, print a message and continue with the next line
        print(f"Error occurred while processing line {line_number}: ", str(e))
        continue

# If in debug mode, print the mutations that have been detected
if debug:
    pprint(detected_mutations)

# If not outputting raw data, print the VCF header
if not rawoutput:
    print("##fileformat=VCFv4.1")
    print("##source=INDELseek")
    print(f"##reference=file://{refseq}")
    print("##INFO=<ID=DP2,Number=2,Type=Integer,Description=\"# alt-foward and alt-reverse reads\">")
    print("##INFO=<ID=QS2,Number=2,Type=Float,Description=\"Mean quality scores of alt-foward and alt-reverse bases\">")
    
    # If a depth file has been provided, include additional header lines
    if depth_bam is not None and depth_bam != "":
        print("##INFO=<ID=AF,Number=1,Type=Float,Description=\"Allele Frequency\">")
        print("##INFO=<ID=RDP,Number=1,Type=Float,Description=\"Mean read depth of REF positions\">")
        print(f"##FILTER=<ID=LowAF,Description=\"AF below {min_af}\">")
    
    print(f"##FILTER=<ID=LowDepth,Description=\"ALT depth below {min_depth}\">")
    print(f"##FILTER=<ID=LowQual,Description=\"Mean quality scores below {quality_threshold} or ALT contains N\">")
    print("\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]))

    # Sort the mutations by their total depth
    mutation_keys = sorted(detected_mutations.keys(), key=lambda k: detected_mutations[k][0] + detected_mutations[k][1], reverse=True)

    # Iterate over each mutation
    for key in mutation_keys:
        # If in debug mode, print the key of the mutation
        if debug:
            print(f"key: {key}")

        # Extract the relevant information from the key
        chrom, start, end, ref, alt = key.split("|")
        
        # Calculate the mean quality score in each direction and overall
        qualscore_forward = mean(detected_mutations[key][2])
        qualscore_reverse = mean(detected_mutations[key][3])
        qualscore_combined = mean(detected_mutations[key][2] + detected_mutations[key][3])
        
        # Calculate the total depth of the mutation
        combined_depth = detected_mutations[key][0] + detected_mutations[key][1]

        # Initialize a list to store filter tags and a string to store the INFO field
        filter_tags = []
        info = ""

        # If the combined quality score is below the threshold or the mutation includes an 'N', append the 'LowQual' filter
        if qualscore_combined < quality_threshold or "N" in alt:
            if skip_lowqual:
                continue
            filter_tags.append("LowQual")
        
        # If the total depth is below the threshold, append the 'LowDepth' filter
        if combined_depth < min_depth:
            if skip_lowdepth:
                continue
            filter_tags.append("LowDepth")

        # Start constructing the INFO field
        info = f"DP2={detected_mutations[key][0]},{detected_mutations[key][1]};QS2={qualscore_forward:.2f},{qualscore_reverse:.2f}"

        # If a depth file has been provided, calculate the allele frequency and mean read depth and add them to the INFO field
        if depth_bam is not None and depth_bam != "":
            region = f"{chrom}:{start}-{end}"
            pos_depth = depth(region)
            af = combined_depth / pos_depth
            if pos_depth == 0:
                raise Exception(f"ERROR: unexpected pos_depth {pos_depth} for {region}")
            if af < min_af:
                if skip_lowaf:
                    continue
                filter_tags.append("LowAF")
            info += f";AF={af:.3f};RDP={pos_depth:.1f}"

        # Print the VCF line for this mutation
        print("\t".join([
            chrom,
            start,
            ".",
            ref,
            alt,
            f"{qualscore_combined:.2f}",
            "PASS" if len(filter_tags) == 0 else ";".join(sorted(filter_tags)),
            info
        ]))

# Get the current time to calculate the runtime of the script
end_time = time.time()

# Calculate and print the runtime of the script
elapsed_time = time.time() - start_time
elapsed_hours = elapsed_time / 60
print(f"\n\nThe script ran in: {elapsed_hours:.2f} minutes")