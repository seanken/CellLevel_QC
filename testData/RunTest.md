## Test Data For Pipeline

We have included a small dataset (3' 10X PMBC single cell data, 100,000 reads) in the outs directory. To test the jar in this directory can run: 

```
java -jar ../Jar/SingleCellQC.jar -d outs -o out.txt -v -s
```

This should produce a file out.txt with the data of interest. The output to consule should look like:

```
Inputs:
Bam: outs/possorted_genome_bam.bam
Cell File: outs/raw_feature_bc_matrix/barcodes.tsv.gz
Output: out.txt
Verbose: false
Cells are gzipped: true
Matrix Directory: outs/raw_feature_bc_matrix
Read in Cell Data
Read in data!
Run time for processing bam: PT1.689S
Total number of alignments: 100000
Process MM Matrix
Not yet implemented
Save file
Run sanity check
Check integer
Percents check out, as do counts!
Compare to metric file
From metric csv, Number of Reads: 100000.0
From counter, total: 100000.0
 
From metric csv, Reads Mapped Confidently to Intergenic Regions: 3.8
From counter, intergenic: 3.796
 
From metric csv, Reads Mapped Confidently to Intronic Regions: 24.2
From counter, intronic: 24.183
 
From metric csv, Reads Mapped Confidently to Exonic Regions: 53.5
From counter, exonic: 53.49
 
From metric csv, Reads Mapped Antisense to Gene: 0.9
From counter, antisense: 0.939
 
From metric csv, Q30 Bases in Barcode: 94.4
From counter, percent_qual_cbc: 94.43581
 
From metric csv, Q30 Bases in UMI: 93.8
From counter, percent_qual_umi: 93.82268
 
Some unit tests
Passed XF test!
Pass splice test!
```
