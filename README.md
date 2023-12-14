# Ecological Genomics Final Project Notebook
UVM PBIO 6800 Ecological Genomics course Fall 2023

### December 8, 2023

- `per_contig_v2` folder has 91,354 files. 15,338 are `_result.Observed.txt` files, moved to `per_contig_results_v2` folder.
- Re-ran all below steps from Dec 6th on new output. Got ~1600 significant contigs for `p.adj <= 0.0005` which is kinda a lot.
- Now need to grep gene IDs and put them into Plant Genie for GO enrichment.

```bash
ANNOT="/netfiles/ecogen/PopulationGenomics/ref_genome/annotation/Pabies1.0-all-cds.gff3.gz"

OUTLIERS=~/finalProject/myresults/signif_contigs_name_only.txt

zcat ${ANNOT} | grep -f ${OUTLIERS} | grep "gene" | uniq | cut -f9 | sed "s/ID=//g" > gene_IDs.txt
```

### December 6, 2023

- `per_contig` folder now has 104,346 files, which when divided by 4 = 26087, as compared to total number of contigs (33678). For each contig, I should in theory have a `.Angsd.arg` file, a `.Angsd.abbababa2` file, a `_result.Observed.txt`, and a `_result.TransRem.txt`. However, I may have the first two for some and not the second two if it couldn’t generate the results files from the R script. Should perhaps move all the `_result.Observed.txt` files to a separate folder.

```bash
# from per_contig directory
mkdir ../per_contig_results
mv *_result.Observed.txt ../per_contig_results
```

- 21834 `_result.Observed.txt` files, 65% of original number of contigs.
- Rename all files to get rid of the irritating colon:

```bash
for f in *:_result.Observed.txt; do mv -- "$f" "${f%:_result.Observed.txt}"_result.Observed.txt; done
```

- The combo we want (H1 = South_RS, H2 = North_RS, H3 = BS, and H4 = NS is on line 5 of every file.

```bash

# put row5 of every file into a new file in a folder

mkdir per_contig_row5 # in results folder

# from per_contig_results folder
for FILE in *_result.Observed.txt; do sed '5q;d' ${FILE} > ../per_contig_row5/${FILE}_row5.txt; done # change row number to get the one you want

# from per_contig_row5 folder
for f in *_result.Observed.txt_row5.txt; do mv -- "$f" "${f%_result.Observed.txt_row5.txt}"_row5.txt; done
```

- If value in column 5 is ≤ 0.05, print file name and then contents of file into one big file.
- Successfully, made file of all contigs with their ABBA-BABA data for the particular combo I wanted with script `combine_files.txt`. For multiple testing correction, should it be the 21834 tests that made it through the initial analysis or the 33678 that includes ones that probably got thrown out for window mismatch (either because window size was bigger than contig or not small enough to break down the contig into sufficient windows for the ABBA-BABA jack-knife. Steve says don’t include ones that got thrown out in multiple testing correction.
- Try with 5000 window size, running now. FYI, even if bam list should only include the `sorted.rmdup.bam` files, the `sorted.rmdup.bam.bai` index files must also be in the same directory, or ANGSD freaks out.

```bash
#!/bin/bash

INPUT=~/finalProject/mydata
OUTPUT=~/finalProject/myresults/per_contig

REF="/netfiles/ecogen/PopulationGenomics/ref_genome/Pabies1.0-genome_reduced.fa"

while read CONTIG

do

ANGSD -doAbbababa2 1 -bam ${INPUT}/RSBS_new_bam.list -r ${CONTIG} -blockSize 5000 -sizeFile ${INPUT}/sizeFile.size -doCounts 1 -out ${OUTPUT}/bam_${CONTIG}.Angsd -anc ${REF} -useLast 0 -minQ 20 -minMapQ 30 -p 1 -nThreads 10

Rscript ~/finalProject/myscripts/estAvgError.R angsdFile=${OUTPUT}/bam_${CONTIG}.Angsd out=${OUTPUT}/${CONTIG}_result sizeFile=${INPUT}/sizeFile.size nameFile=${INPUT}/popNames.names

done < ${INPUT}/all_contigs.txt # used a test_contigs.txt folder with only 10 contigs for early troubleshooting

# [E::idx_find_and_load] Could not retrieve index file for '/netfiles/ecogen/groupProjects/sprucedUp/RS_mapped_bam/2019_25_N.sorted.rmdup.bam'
# [main_samview] random alignment retrieval only works for indexed BAM or CRAM files. (file: '/netfiles/ecogen/groupProjects/sprucedUp/RS_mapped_bam/2019_25_N.sorted.rmdup.bam' )
```

### December 5, 2023

- Re-run overall ABBA-BABA with new bam files mapped to Norway Spruce genome and lump North_RS and Mid_RS together.

```bash

ls "/netfiles/ecogen/groupProjects/sprucedUp/RS_mapped_bam"/*sorted.rmdup.bam > RSBS_new_bam.list

# removed Steve's 444 test population manually in Atom
```

- Figure out next step for per-contig analysis. Need to write code to go through all the files and pull out ones with a significant p-value for the comparison where (from left to right) H1 = South_RS, H2 = Mid_North_RS, H3 = BS, and H4 = NS.
- Re-run per-contig analysis with above parameters.
- Figure out how to run q-value correction for multiple testing in R from a list of p-values.

### December 4, 2023

- Conducted ABBA-BABA window analysis with window sizes 10k, 5k, and 1k bp. Got identical outputs in R, except for numerical variation in `Z` and `nBlocks`. All are still significant.
- Found these tutorials ([1](https://github.com/simonhmartin/tutorials/blob/master/ABBA_BABA_windows/README.md)) ([2](https://evomics.org/learning/population-and-speciation-genomics/2018-population-and-speciation-genomics/abba-baba-statistics/)) on ABBA-BABA sliding windows that use Python for whole-genome ABBA-BABA and give more of an output like I’d expect.
- From Norway spruce ref genome FASTA file, made list of contig names with colon appended after each.
    - `grep “MA” or “>” | sed ‘s/>/g’ | sort | uniq | head`
- Make a while read loop, to run ABBA-BABA step using the regular ANGSD `-rf [region file]` flag, to do the ABBA-BABA on each window one at a time. Will also need to set block size within each contig to something reasonable.

```bash
#!/bin/bash

INPUT=~/finalProject/mydata
OUTPUT=~/finalProject/myresults/per_contig

REF="/netfiles/ecogen/PopulationGenomics/ref_genome/Pabies1.0-genome_reduced.fa"

while read CONTIG

do

ANGSD -doAbbababa2 1 -bam ${INPUT}/RSBS_bam.list -r ${CONTIG} -blockSize 20 -sizeFile ${INPUT}/sizeFile.size -doCounts 1 -out ${OUTPUT}/bam_${CONTIG}.Angsd -anc ${REF} -useLast 0 -minQ 20 -minMapQ 30 -p 1 -nThreads 5

Rscript ~/finalProject/myscripts/estAvgError.R angsdFile=${OUTPUT}/bam_${CONTIG}.Angsd out=${OUTPUT}/${CONTIG}_result sizeFile=${INPUT}/sizeFile.size nameFile=${INPUT}/popNames.names

done < ${INPUT}/all_contigs.txt

# Error in seq.default(1, lenList, numComb) : wrong sign in 'by' argument
# Calls: seq -> seq.default
# Execution halted

# Error in `colnames<-`(`*tmp*`, value = colnames(outDataTotal)) : 
#   attempt to set 'colnames' on an object with less than two dimensions
# Execution halted
```

- Above code is in running and producing below error, but I think it’s getting tripped up on the error correction step. Output files appear ok. 33678 contigs total. Chose block size of 20 based on ball-parking: ANGSD defaults were using 5 Mbp blocks on the ~1.3 Gbp reduced reference genome = ~275 blocks for the whole genome. So if each contig is ~5 kbp, blocks of size 20 bp would give us ~250 blocks.
    - Steve says this is too small because the genome-wide SNP rate is 3 per 1000 bp. The `seq.default` error may be because some of my blocks have no SNPs in them. Increased to 1500 and restarted run. This creates a new error on some, but not all contigs. Tentatively, contigs with errors do not seem to be writing to the results folder, which is good.
- I have drastically increased the number of tests by doing this, and I wonder if I’ll need to correct the p-values myself.
    - `qqman` q-value correction R package, false discovery based correction, less severe than Bonferonni.
- Next step: Pick a specific combination of populations, pull that row out of each file, only keep contigs that are significant for that combo. Should I lump North_RS and Mid_RS together in that case?
- Outlier contigs from genetic PCA are in Gwen’s GitHub.
- Run ABBA-BABA on new bam files for Norway spruce mapping.

### November 29, 2023

- Removed `errFile=errorList.error` in below R command. Made popNames.name file:
    - South_RS
    - Mid_RS
    - North_RS
    - BS
    - NS

```jsx
#!/bin/bash

DIR=~/finalProject/mydata

Rscript ~/finalProject/myscripts/estAvgError.R angsdFile=${DIR}/bam.Angsd.abbababa2 out=${DIR}/result sizeFile=${DIR}/sizeFile.size nameFile=${DIR}/popNames.names
```

- Downloaded `estAvgError.R` from [https://github.com/ANGSD/angsd/blob/master/R/estAvgError.R](https://github.com/ANGSD/angsd/blob/master/R/estAvgError.R)
    - Csenge needed to use her sudo powers to install the `pacma` R library.
- Outputted files:
    - `result.Observed.txt` = D-statistic calculated WITHOUT Error Correction and WITHOUT Ancient Transition removal
    - `result.TransRem.txt` = D-statistic calculated WITHOUT Error Correction and WITH Ancient Transition removal (probably not important)
- All comparisons very significant, but given an alternative, at no point did it find introgression between black spruce and the south. It’s all significant because the hybridization was recent and the lineages are closely related, diverged recently, and are not particularly reproductively isolated. Interpreting magnitude of the D-statistic may be meaningful.
- Now try with window size of 5k and 10k, using `blockSize` field in ANGSD, combined ANGSD and R script into one step.

```jsx
#!/bin/bash

DIR=~/finalProject/mydata

REF="/netfiles/ecogen/PopulationGenomics/ref_genome/Pabies1.0-genome_reduced.fa"

ANGSD -doAbbababa2 1 -bam ${DIR}/RSBS_bam.list -blockSize 10000 -sizeFile ${DIR}/sizeFile.size -doCounts 1 -out ${DIR}/bam_window.Angsd -anc ${REF} -useLast 0 -minQ 20 -minMapQ 30 -p 1

Rscript ~/finalProject/myscripts/estAvgError.R angsdFile=${DIR}/bam_window.Angsd.abbababa2 out=${DIR}/window_result sizeFile=${DIR}/sizeFile.size nameFile=${DIR}/popNames.names
```

- Moved done things into shared folder: `/netfiles/ecogen/groupProjects/sprucedUp`

### November 27, 2023

- ABBA-BABA multipop test: [http://popgen.dk/angsd/index.php/Abbababa2](http://popgen.dk/angsd/index.php/Abbababa2)
    - Copied RSBS_bam.list to my directory
    - Going to use indexed P. abies ref genome: `/netfiles/ecogen/PopulationGenomics/ref_genome/Pabies1.0-genome_reduced.fa.fai`
    - Skipped error correction step. This is only for ancient DNA to correct for transition errors where cytosine switches to thymine due to some degenerative process over time.
- sizeFile option specifies number of individuals in each population (more than 4 populations can be defined). If not provided, it is assumed that each population has only one individual. We have 95 individuals across 12 populations, 8 per red spruce pop. except 2100 which only has 7. For black spruce pops, we have: MN1 = 4, MN2 = 5, WISC = 9.
- Only want 5 "populations": south RS (2019, 2020), middle RS (2021), north RS, BS, and outgroup, so I need to pool them like that. Then ABBA-BABA will test every 4 population combination of those.
    - south RS: 16
    - middle RS: 8
    - north RS: 71
    - BS: 18
    - outgroup: 1

```jsx
#!/bin/bash

DIR=~/finalProject/mydata

REF="/netfiles/ecogen/PopulationGenomics/ref_genome/Pabies1.0-genome_reduced.fa"

ANGSD -doAbbababa2 1 -bam ${DIR}/RSBS_bam.list -sizeFile ${DIR}/sizeFile.size -doCounts 1 -out ${DIR}/bam.Angsd -anc ${REF} -useLast 0 -minQ 20 -minMapQ 30 -p 1
```
