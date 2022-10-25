# CellLevel_QC

This repository contains java code for a simple command to get cell level mapping information (% reads intronic, intergenic, multimapped, etc) from the output in CellRanger (tested with CellRanger v6, should work with more recent versions as well though no gurantees).

## How to use

To get the cell level QC information use Jar/SingleCellQC.jar. The script requires java (was compiled with version 1.8). 

The most basic way to run, if input.bam is the input bam, cells.txt a list of cell barcodes (plain text format, for gzipped see below, cell names should have the -1 suffix), out.txt is the file you want your results to be writen to, then can run

'''
java -jar /path/to/repo/Jar/SingleCellQC.jar -i input.bam -c cells.txt -o out.txt
'''

If you'd like a more verbose version (produces and output every 1 million reads) you can add the -v flag:

'''
java -jar /path/to/repo/Jar/SingleCellQC.jar -i input.bam -c cells.txt -o out.txt -v
'''

If your cell file is gzipped (so cells.txt.gz) as the list produced by CellRanger is need the -z flag:

'''
java -jar /path/to/repo/Jar/SingleCellQC.jar -i input.bam -c cells.txt.gz -o out.txt -z
'''
