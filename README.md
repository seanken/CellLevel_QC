# CellLevel_QC

This repository contains java code for a simple command to get cell level mapping information (% reads intronic, intergenic, multimapped, etc) from the output in CellRanger (tested with CellRanger v6, should work with more recent versions as well though no gurantees).

## How to use

To get the cell level QC information use Jar/SingleCellQC.jar. The script requires java (was compiled with version 1.8). 

The most basic way to run, if input.bam is the input bam, cells.txt a list of cell barcodes (plain text format, for gzipped see below, cell names should have the -1 suffix), out.txt is the file you want your results to be writen to, then can run

```
java -jar /path/to/repo/Jar/SingleCellQC.jar -i input.bam -c cells.txt -o out.txt
```

If you'd like a more verbose version (produces and output every 1 million reads) you can add the -v flag:

```
java -jar /path/to/repo/Jar/SingleCellQC.jar -i input.bam -c cells.txt -o out.txt -v
```

If your cell file is gzipped (so cells.txt.gz) as the list produced by CellRanger is need the -z flag:

```
java -jar /path/to/repo/Jar/SingleCellQC.jar -i input.bam -c cells.txt.gz -o out.txt -z
```

Alternativelly, can feed in a 10X out directory, which will find the bam file, cell file (in the raw directory), and matrix (not used atm):

```
java -jar /path/to/repo/Jar/SingleCellQC.jar -d /path/to/cellranger/output/outs
```

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

`polyA:` Number of reads with any number of bases trimmed due to poly-A (does not count unmapped reads).

`TSO:` Number of reads with any number of bases trimmed due to TSO (does not count unmapped reads).

`spliced: ` Number of reads with a splice even in them (only counts uniquelly mapped reads)

`percent_qual_cbc: ` Average percent of bases (rounded to the nearest percent) in the CBC over all reads with quality score >30.

`percent_qual_umi: ` Average percent of bases (rounded to the nearest percent) in the UMI over all reads with quality score >30.

`total:` The total number of reads in this cell.

Will likely add more/try to extend beyond CellRanger but that might not be for some time.

## Build

The jar file was built with gradle. To rerun gradle change directory to the gradle directory and run:

```
./gradlew build
```

The new jar will be in app/build/libs.

## Code

The src code is in the src directory. Will work on making cleaner code if I get the chance.
