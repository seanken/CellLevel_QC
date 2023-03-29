package singlecellqc;
import java.util.*;
import java.lang.*;
import org.apache.commons.cli.*;


//
//This is the class that extracts the input arguments
//
public class CellLevel_QC
{
    public static void main(String[] args)
    {
        CommandLine cmd=processInputs(args);
        GetQC(cmd);

    }

    //processes the inputs so can be used to run command
    public static CommandLine processInputs(String[] args)
    {
    
        //Sets up options to use
        Options options = new Options();
        Option inputBam = new Option("i", "input", true, "input BAM file from CellRanger (overriden by -d)");
        options.addOption(inputBam);

        Option inputCell = new Option("c", "cells", true, "input cell names (overriden by -d)");
        options.addOption(inputCell);

        Option quantUsed = new Option("q", "quantused", true, "quantification method used (CellRanger or STARSolo). If STARSolo is used preprocessing is required (see github). Default is CellRanger.");
        options.addOption(quantUsed);

        Option inputMat = new Option("m", "matrix", true, "input UMI count matrix in MM format as in CellRanger (overridden by -d, not currently used)");
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

        Option allReads = new Option("a", "all", false, "includes multimappers in all metrics");
        options.addOption(allReads);

        Option gzipped = new Option("z", "gzipped", false, "if cells are gzipped (set to true if -d is given)");
        options.addOption(gzipped);

        Option testing = new Option("t", "test", false, "if want to run in testing mode (only use first 10,000,000 reads if true)");
        options.addOption(testing);

        Option checking = new Option("s", "sanityCheck", false, "runs a few sanity check on results (only works with -d argument)");
        options.addOption(checking);

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



    //Based on the inputs runs the QC commmands
    public static void GetQC(CommandLine cmd)
    {
        String inputBamPath=null; //the input bam from CellRanger or STARSolo
        String inputCellPath=null; //The path the the file with cell names in it
        String inputMatPath=null; //The path to the MM matrix
        String inputGTFPath=null; //Path to the gtf used by CellRanger or STARSolo (optional)
        String indirPath=null; //Input CellRanger directory if passed as argument
        String quantUsed="CellRanger"; //Tells us if CellRanger or STARSolo

        
        //Figure out bam/call path for later use
        if(!cmd.hasOption("d"))
        {
            if(!cmd.hasOption("i") | !cmd.hasOption("c"))
            {
                print("Need either the -d argument or the -i and -c arguments");
                System.exit(1);
            }
            inputBamPath = cmd.getOptionValue("input");
            inputCellPath = cmd.getOptionValue("cells");
            if(cmd.hasOption("m")){inputMatPath = cmd.getOptionValue("matrix");};
        }
        else
        {
            indirPath=cmd.getOptionValue("d");
            inputBamPath=indirPath+"/possorted_genome_bam.bam";
            inputCellPath=indirPath+"/raw_feature_bc_matrix/barcodes.tsv.gz";
            inputMatPath=indirPath+"/raw_feature_bc_matrix";
        }        


        //Checks what quantification method is used
        if(cmd.hasOption("q"))
        {
            quantUsed=cmd.getOptionValue("quantused");
            
            if(!quantUsed.equals("STARSolo") & !quantUsed.equals("CellRanger"))
            {
                print("Quantification method used (-q option) must be CellRanger or STARSolo");
                return;
            }
        }

        if(quantUsed=="STARSolo" & cmd.hasOption("d"))
        {
            print("Option -d can not be used with STARSolo");
            return;
        }





        String outputPath = cmd.getOptionValue("output"); //Text file output QC metrics to

        //boolean options
        boolean verboseVal=cmd.hasOption("v");
        boolean gzipCells=cmd.hasOption("z");
        boolean useMulti=cmd.hasOption("a");

        if(cmd.hasOption("d"))
        {
            gzipCells=true;
        }
        boolean testingVal=cmd.hasOption("t");


        //Prints out settings to screen
        print("Inputs:");
        print("Bam: "+inputBamPath);
        print("Cell File: "+inputCellPath);
        print("Output: "+outputPath);
        print("Verbose: "+String.valueOf(verboseVal));
        print("Use multimappers: "+String.valueOf(useMulti));
        print("Cells are gzipped: "+String.valueOf(gzipCells));
        print("Quantification method: "+quantUsed);
        if(inputMatPath==null)
        {
            print("Matrix Directory: none");
        }else{
            print("Matrix Directory: "+inputMatPath);
        }

        //Creates ReadCounter object that will be used for extracting cell level QC with given arguments
        ReadCounter counter=new ReadCounter(inputBamPath,inputCellPath,outputPath,gzipCells,quantUsed,useMulti);
        
        //If gtf is given processes it
        if(cmd.hasOption("g"))
        {
            inputGTFPath=cmd.getOptionValue("gtf");
            print("Process GTF");
            counter.ProcessGTF(inputGTFPath);
        }


        counter.ReadBam(verboseVal,testingVal); //Processes the bam to get QC metrics
        
        //if(inputMatPath!=null)
        //{
        //    print("Process MM Matrix");
        //    counter.ProcessMatrix(inputMatPath);
        //}
        //else{
        //    print("No MM Matrix given, skipping processing");
        //}
        counter.SaveQC(); //save the QC results to table

        //For testing
        if(cmd.hasOption("s") & cmd.hasOption("d"))
        {
            String metricsComb=indirPath+"/metrics_summary.csv";
            print("Run sanity check");
            TestCounter testCount=new TestCounter();
            print("Check integer");
            testCount.checkInteger(counter);
            print("Compare to metric file");
            testCount.compareToMetricCSV(counter,metricsComb);
            print("Some unit tests");
            testCount.checkXFParsing(counter);
            testCount.checkSplice(counter);
            //to implement
            testCount.checkRegionMappingTo(counter);
            //testCount.checkTrim(counter);
            print("To be added: UTR unit tests");

        }

    }


    public static void print(String output)
    {
        System.out.println(output);
    }

}
