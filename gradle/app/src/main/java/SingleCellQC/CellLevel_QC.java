package SingleCellQC;
import java.util.*;
import java.lang.*;
import org.apache.commons.cli.*;



public class CellLevel_QC
{
    public static void main(String[] args)
    {
        Options options = new Options();
        Option inputBam = new Option("i", "input", true, "input BAM file from CellRanger");
        inputBam.setRequired(true);
        options.addOption(inputBam);

        Option inputCell = new Option("c", "cells", true, "input cell names");
        inputCell.setRequired(true);
        options.addOption(inputCell);

        Option output = new Option("o", "output", true, "output file path");
        output.setRequired(true);
        options.addOption(output);

        //Option inputDS = new Option("p", "prop", true, "proportion to downsample to, 1 by default (no DS)");
        //options.addOption(inputDS);

        Option verbose = new Option("v", "verbose", false, "to make verbose");
        options.addOption(verbose);

        Option gzipped = new Option("z", "gzipped", false, "if cells are gzipped");
        options.addOption(gzipped);

        Option testing = new Option("t", "test", false, "if want to run in testing mode (only use first 10,000,000 reads)");
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

        String inputBamPath = cmd.getOptionValue("input");
        String inputCellPath = cmd.getOptionValue("cells");
        String outputPath = cmd.getOptionValue("output");

        //double DS=1;
        //if(cmd.hasOption("p"))
        //{
        //    try{
        //        DS=Double.parseDouble(cmd.getOptionValue("prop"));
        //    }catch(Exception e)
        //    {
        //        print("Not valid prop argument");
        //        System.exit(0);
        //    }
        //}

        boolean verboseVal=cmd.hasOption("v");
        boolean gzipCells=cmd.hasOption("z");
        boolean testingVal=cmd.hasOption("t");

        print("Inputs:");
        print("Bam: "+inputBamPath);
        print("Cell File: "+inputCellPath);
        print("Output: "+outputPath);
        print("Verbose: "+String.valueOf(verboseVal));
        ReadCounter counter=new ReadCounter(inputBamPath,inputCellPath,outputPath,gzipCells);
        counter.ReadBam(verboseVal,testingVal);
        counter.SaveQC();

    }


    public static void print(String output)
    {
        System.out.println(output);
    }

}
