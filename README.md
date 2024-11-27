
# INDELseek Algorithm Reimplementation and Optimization

This repository contains the implementation and optimization of the INDELseek algorithm, an open-source tool designed for detecting complex insertions and deletions (indels) in Next-Generation Sequencing (NGS) data, as described the ["INDELseek: detection of complex insertions and deletions from next-generation sequencing data" by Au C.H. et. al](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-016-3449-9). Originally implemented in Perl, INDELseek is a crucial tool for identifying genetic variations that may play significant roles in diseases, including cancer. Our project focuses on re-implementing the algorithm in Python and C++ to enhance performance and usability, as well as comparing it to existing tools like GATK.

This analysis was performed as part of the final project for the "Algorithms in Molecular Biology" graduate course of the MSc Data Science & Information Technologies Master's programme (Bioinformatics - Biomedical Data Science Specialization) of the Department of Informatics and Telecommunications department of the National and Kapodistrian University of Athens (NKUA), under the supervision of professors Theodore Dalamagas and George Georgakilas, in the academic year 2022-2023.

---

## Contributors
- [Konstantinos Giatras](https://github.com/GiatrasKon)
- [Alexandros Kouris](https://github.com/Alx-Kouris)

---

## Workflow and Tools

### Key Components
- **Original INDELseek Algorithm**: Perl script enhanced with execution time tracking and memory logging.
- **Python Implementation**: Uses dictionaries for caching and simplifies debugging with modular design.
- **C++ Implementation**: Incorporates hash maps for caching, custom structs for memory efficiency, and a `Makefile` for easy compilation.
- **Memory Usage Visualization**: `graph.sh` generates graphs for RAM usage over time during execution.

### Tools and Datasets
- **Datasets**: NGS datasets including NA12878 (1.5 billion reads) and HG00148.chr11 (4 million reads) were used.
- **Computing Environment**: Testing was performed locally and on a Hypatia VM (7 cores, 60GB RAM).
- **Dependencies**: Scripts rely on `samtools` for sequence processing and alignment.

---

## Results

- **Execution Time**: The C++ implementation demonstrated significant speed improvements (40 minutes for ~4M reads) compared to Python (5 hours 42 minutes) and Perl (5 hours 49 minutes).
- **Memory Usage**: Resource footprint was consistently low (~600MB RAM for a 4M reads dataset).
- **Insights**: Compiled languages like C++ show clear performance advantages, albeit at increased coding complexity.

---

## Installation

### Cloning the Repository

```sh
git clone https://github.com/GiatrasKon/INDELseek-Optimized-Implementations.git
```

### Required Tools

- `samtools`
    - Processes BAM/SAM files for read alignments.
    - Fetches reference sequences for specific regions using `faidx`.
    - Computes depth of coverage using `depth`.
    - Installation:
```sh
sudo apt install samtools
```
- `gnuplot`
    - Used in `graph.sh` to generate memory usage graphs for the running process.
    - Installation:
```sh
sudo apt install gnuplot
```
- `perl`
    - The original INDELseek implementation is written in Perl. If you plan to compare outputs, ensure Perl is installed with the required modules like `Getopt::Long`.
    - Installation:
```sh
sudo apt install perl
```


### Dependencies

**Python**:

Ensure you have Python 3.x installed. Use the following commands to install required packages:

```sh
pip install argparse psutil subprocess datetime time collections pprint
```

**C++ and Perl**:

- Ensure `gcc` (C++ compiler) and `make` are installed.
- Install `samtools` (v1.3 or higher).
- Verify Perl is installed with required modules like `Getopt::Long`.
- **Compile the C++ Implementation**: Navigate to the scripts directory and run:

```sh
make
```

### Datasets

1. **NA12878 BAM File**:
    - This is a 113GB BAM file containing approximately 1.5 billion reads. It is aligned to the hg38 reference genome.
    - Used to attempt reproduction of the original INDELseek results.
    - Extracted subsets (chromosomes 20 and 15) for faster testing and analysis.
    - How to Obtain:
        - Visit [Illumina BaseSpace](https://login.illumina.com/platform-services-manager/?rURL=https://basespace.illumina.com&clientId=basespace&clientVars=aHR0cHM6Ly9iYXNlc3BhY2UuaWxsdW1pbmEuY29tL2Rhc2hib2FyZA&redirectMethod=GET#/) (requires an account).
        - Look for the NA12878 sample in publicly available datasets or the [Genome in a Bottle (GIAB) project](https://www.nist.gov/programs-projects/genome-bottle).
2. **HG00148.chr11 BAM File**:
    - A BAM file containing approximately 4 million reads focused on chromosome 11.
    - Used for benchmarking the Python and C++ implementations.
    - Provides a smaller, manageable dataset for performance comparisons.
    - How to Obtain:
        - Visit the [1000 Genomes Project website](https://www.internationalgenome.org/).
        - Download the HG00148 dataset.
3. **Reference Genome (hg38)**:
    - A FASTA file of the human genome (hg38 version) used for sequence alignment and comparison.
    - Required for samtools operations like faidx and depth.
    - Essential for running all INDELseek implementations.
    - How to Obtain:
        - [Broad Institute Public Repository](https://www.broadinstitute.org/datasets).
        - [UCSC Genome Browser](https://genome.ucsc.edu/).

How to Download and Prepare the Datasets:
1. NA12878 and HG00148:
        - BAM files can be downloaded from their respective sources using web interfaces or wget/curl commands (if direct URLs are provided).
2. Reference Genome:
        - Download the hg38.fa file and index it using samtools:
```sh
samtools faidx hg38.fa
```

Notes for Users:
- Ensure sufficient storage space for large datasets.
- Validate downloaded files (e.g., by comparing checksums).
- Ensure compatibility between the reference genome version and BAM file alignments.
- Follow ethical guidelines and licensing agreements when using public datasets.

---

## Usage

Python:

```sh
python indelseek.py --ref_genome path/to/reference.fa --input_bam path/to/input.bam
```

C++:

```sh
./indelseek --reference path/to/reference.fa --input_bam path/to/input.bam --output_vcf output.vcf
```
---

## Documentation

Refer to the `documents` directory for the project presentation, report, and the original authors' paper.

---