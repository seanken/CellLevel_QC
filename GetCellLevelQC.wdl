version 1.0

workflow GetCellLevelQC {
    input {
        File? in_dir
        File bam_file=in_dir+"/possorted_genome_bam.bam" ##A BAM file containing aligned reads from single-cell RNA-seq experiment, from the CellRanger pipeline
        File cells_file=in_dir+"/filtered_feature_bc_matrix/barcodes.tsv.gz" ##A file containing cell barcodes to be analyzed, expected to be gzipped
        String output_prefix = "cell_qc"
    }

    call RunQC {
        input:
            bam_file = bam_file,
            cells_file = cells_file,
            output_prefix = output_prefix
    }

    output {
        File qc_output = RunQC.qc_results
    }
}

task RunQC {
    input {
        File bam_file
        File cells_file
        String output_prefix
        
        # Runtime parameters
        Int memory_gb = 30
        Int disk_space_gb = 100
        Int cpu = 1
        String docker = "amazoncorretto:17"
    }

    String output_filename = "~{output_prefix}.txt"
    String jar_filename = "qc_tool.jar"

    command <<<
        # Download JAR file from GitHub
        wget https://github.com/seanken/CellLevel_QC/raw/refs/heads/main/Jar/SingleCellQC.jar
        
        # Run the QC tool
        java -jar SingleCellQC.jar \
            -i ~{bam_file} \
            -c ~{cells_file} \
            -o ~{output_filename} \
            -z
    >>>

    output {
        File qc_results = output_filename
    }

    runtime {
        docker: docker
        memory: "~{memory_gb} GB"
        disks: "local-disk ~{disk_space_gb} HDD"
        cpu: cpu
    }
}
