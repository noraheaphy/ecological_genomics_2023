# Population Genomics Lab Notebook

## Author: Nora Heaphy

### Affiliation: University of Vermont, Dept. of Plant Biology

### E-mail contact: [nora.heaphy\@uvm.edu](mailto:nora.heaphy@uvm.edu){.email}

### Start Date: Sept. 11, 2023

### End Date:

### Project Descriptions: This notebook will document my workflow on the bioinformatics of the Population Genomics section of the Fall 2023 Ecological Genomics course.

# Table of Contents:

-   [Entry 1: 2023-09-11](#id-section1)
-   [Entry 2: 2023-09-13](#id-section2)
-   [Entry 3: 2023-09-18](#id-section3)

------------------------------------------------------------------------

<div id='id-section1'/>

### Entry 1: 2023-09-11.

-   We reviewed the red spruce study system and the exome capture data.
-   We discussed the structure of fastq files (DNA sequences, Q-scores).
-   Using the program FastQC, we analyzed the quality of the sequencing runs for one file.

------------------------------------------------------------------------

<div id='id-section2'/>

### Entry 2: 2023-09-13.

-   After discussing FastQC results, we saw good quality sequence data for most of the read length.
-   The initial 5 bp or so had more variable base frequencies, and the very end of the reads had slightly lower Q-scores.
-   Based on this, we set up an analysis to trim the reads using the 'fastp' program.
-   We ran the bash script 'fastp.sh' for this.
-   We looked at the html files produced by 'fastp' and compared pre- and post-trimming. Things looked good!
-   We ended the day be setting up our read mapping of the trimmed and cleaned reads using 'bwa'.

------------------------------------------------------------------------

<div id='id-section3'/>

### Entry 3: 2023-09-18.
