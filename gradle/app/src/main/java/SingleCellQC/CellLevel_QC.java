package SingleCellQC;
import java.util.*;
import java.lang.*;
import org.apache.commons.cli.*;



public class CellLevel_QC
{
    public static void main(String[] args)
    {
        CommandLine cmd=processInputs(args);
        GetQC(cmd);

    }

    public static CommandLine processInputs(String[] args)
    {
    
        Options options = new Options();
        Option inputBam = new Option("i", "input", true, "input BAM file from CellRanger (overriden by -d)");
        //inputBam.setRequired(true);
        options.addOption(inputBam);

        Option inputCell = new Option("c", "cells", true, "input cell names (overriden by -d)");
        //inputCell.setRequired(true);
        options.addOption(inputCell);

        Option inputMat = new Option("m", "matrix", true, "input UMI count matrix in MM format as in CellRanger (overridden by -d)");
        options.addOption(inputMat);

        Option inputGTF = new Option("g", "gtf", true, "gtf used for CellRanger reference (required for UTR info)");
        options.addOption(inputGTF);

        Option inputDir = new Option("d", "inputDir", true, "input CellRanger outs directory (overrides --input and --cells)");
        options.addOption(inputDir);

        Option output = new Option("o", "output", true, "output file path");
        output.setRequired(true);
        options.addOption(output);


        Option verbose = new Option("v", "verbose", false, "to make verbose");
        options.addOption(verbose);

        Option gzipped = new Option("z", "gzipped", false, "if cells are gzipped (set to true if -d is given)");
        options.addOption(gzipped);

        Option testing = new Option("t", "test", false, "if want to run in testing mode (only use first 10,000,000 reads if true)");
        options.addOption(testing);

        CommandLineParser parser = new DefaultParser();
        HelpFormatter formatter = new HelpFormatter();
        CommandLine cmd=null;

        try {
            cmd = parser.parse(options, args);
        } catch (ParseException e) {
            System.out.println(e.getMessage());
            formatter.printHelp("celllevel_qc", options);

            System.exit(1);
        }
        return(cmd);
    }


    public static void GetQC(CommandLine cmd)
    {
        String inputBamPath=null;
        String inputCellPath=null;
        String inputMatPath=null;
        String inputGTFPath=null;

        if(!cmd.hasOption("d"))
        {
            if(!cmd.hasOption("i") | !cmd.hasOption("c"))
            {
                print("Need either the -d argument or the -i and -c arguments");
                //formatter.printHelp("celllevel_qc", options);
                System.exit(1);
            }
            inputBamPath = cmd.getOptionValue("input");
            inputCellPath = cmd.getOptionValue("cells");
            if(cmd.hasOption("m")){inputMatPath = cmd.getOptionValue("matrix");};
        }
        else
        {
            String indirPath=cmd.getOptionValue("d");
            inputBamPath=indirPath+"/possorted_genome_bam.bam";
            inputCellPath=indirPath+"/raw_feature_bc_matrix/barcodes.tsv.gz";
            inputMatPath=indirPath+"/raw_feature_bc_matrix";
        }        


        String outputPath = cmd.getOptionValue("output");

        boolean verboseVal=cmd.hasOption("v");
        boolean gzipCells=cmd.hasOption("z");
        if(cmd.hasOption("d"))
        {
            gzipCells=true;
        }
        boolean testingVal=cmd.hasOption("t");

        print("Inputs:");
        print("Bam: "+inputBamPath);
        print("Cell File: "+inputCellPath);
        print("Output: "+outputPath);
        print("Verbose: "+String.valueOf(verboseVal));
        print("Cells are gzipped: "+String.valueOf(gzipCells));
        if(inputMatPath==null)
        {
            print("Matrix Directory: none");
        }else{
            print("Matrix Directory: "+inputMatPath);
        }


        ReadCounter counter=new ReadCounter(inputBamPath,inputCellPath,outputPath,gzipCells);
        if(inputGTFPath!=null)
        {
            print("Process GTF");
            counter.ProcessGTF(inputGTFPath);
        }
        counter.ReadBam(verboseVal,testingVal);
        if(inputMatPath!=null)
        {
            print("Process MM Matrix");
            counter.ProcessMatrix(inputMatPath);
        }
        else{
            print("No MM Matrix given, skipping processing");
        }
        counter.SaveQC();

    }


    public static void print(String output)
    {
        System.out.println(output);
    }

}
