# CellLevel_QC: Extracting cell level QC metrics from CellRanger barcoded bams

This repository contains java code for a simple command to get cell level mapping information (% reads intronic, intergenic, multimapped, etc) from the output in CellRanger (tested with CellRanger v6 using the include introns flag, should work with more recent versions as well though no gurantees. Should work on spaceranger output as well, and both single nuc and single cell data).

## How to use

To get the cell level QC information use Jar/SingleCellQC.jar. The script requires java (was compiled with version 1.8). 

The most basic way to run, given an outs directory produced by CellRanger (tested with v6), if you want the output to go to out.txt, then you can run

```
java -jar /path/to/repo/Jar/SingleCellQC.jar -d /path/to/cellranger/output/outs -o out.txt
```

If one would like to count the number of reads overlapping a UTR (counts reads overlapping both the 3' and 5' end) need to pass a gtf with UTRs in it (can have UTR, three_prime_utr, or fivee_prime_utr):

```
java -jar /path/to/repo/Jar/SingleCellQC.jar -d /path/to/cellranger/output/outs -o out.txt -g genes.gtf
```

For a more verbose version (printing out a message every million alignments) can add the -v flag. Can also use the -t flag to only run on the first 10 million reads (for testing purposes). Finally there is also the -s flag that runs some basic sanity checks (checking things that should be integers are, things that should be percentages are, and comparing the results to the sample level (as opposed to cell level) metrics file output by 10X).

If you do not have the outs directory or only want to run on a subset of cells, etc, the most basic way to run, if input.bam is the input bam, cells.txt a list of cell barcodes (plain text format, for gzipped see below, cell names should have the -1 suffix), out.txt is the file you want your results to be writen to, then can run

```
java -jar /path/to/repo/Jar/SingleCellQC.jar -i input.bam -c cells.txt -o out.txt
```

If your cell file is gzipped (so cells.txt.gz) as the list produced by CellRanger is need the -z flag:

```
java -jar /path/to/repo/Jar/SingleCellQC.jar -i input.bam -c cells.txt.gz -o out.txt -z
```

Note -o is required, as is either -d or -c and -i.

A more detailed description will be added, as well examples and tests, assuming I ever have the time.

## Output

The output is a matrix with one row per cell barcode in the cell list given. The columns:

`CBC:` Cell barcode

`antisense:` Number of reads that are antisense to some gene (excludes one that are sense in relation to some gene).

`intergenic:` Number of reads that are intergenic.

`intronic:` Number of reads that are intronic.

`exonic:` Number of reads that are exonic.

`multi:` Number of reads that are multimappers (excluded from exonic/intronic/intergenic numbers).

`unmapped:` Number of reads that are unmapped.

`highConf:` Number of reads that map to transcriptome with high conifdence.

`polyA:` Number of reads with any number of bases trimmed due to poly-A.

`TSO:` Number of reads with any number of bases trimmed due to TSO.

`spliced:` Number of reads with a splice even in them (only counts uniquelly mapped reads)

`percent_qual_cbc:` Average percent of bases (rounded to the nearest percent) in the CBC over all reads with quality score >30.

`percent_qual_umi:` Average percent of bases (rounded to the nearest percent) in the UMI over all reads with quality score >30.

`nUMI:` Number of UMI in this cell.

`UTR:` Number of reads overlapping a known UTR of the gene the are assigned to. Set to 0 if -g is not given.

`total:` The total number of reads in this cell.

Will likely add more/try to extend beyond CellRanger but that might not be for some time.

## Arguments

`-o,--output:` The output file name (required)

`-d,--inputDir:` Path to the outs directory produced by CellRanger. Either need this or -c and -i.

`-c,--cells:` Tab seperated list of cells, one per line, with each cell name ending in -1. Overridden by -d. 

`-z,--gzipped:` Used if the cell list passed to -c is zipped. Set automatically if -d is given.

`-i,--input:` A barcoded bam file produced by CellRanger. Overridden by -d.
 
 `-g,--gtf:` A gtf file, preferably one matching the CellRanger reference. Used to figure out which reads overlap a UTR. If not given UTR is set to 0 for all cells.

 `-v,--verbose:` Included for a more verbose output (print a message every million lines of bam file).

 `-t,--test:` Used for testing purposes, ends the program after reading in 10 million alignments from the bam.

 `-s,--sanityCheck:` Prints some basic sanity checks of the analysis. In particular, looks to make sure integers are integers, percents are percents, and compares the aggregate results of the QC metrics produced to the sample level results returned by CellRanger. Meant to be used with the -d flag.

 `-m,--matrix:` Matrix input, not currently used though likely will be some day.


## Build

The jar file was built with gradle. To rerun gradle change directory to the gradle directory and run:

```
./gradlew build
```

The new jar will be in app/build/libs.

## Code

The src code is in the src directory. Will work on making cleaner code if I get the chance.
