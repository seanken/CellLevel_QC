package SingleCellQC;
import java.util.*;
import java.lang.*;
import java.io.*;

public class TestCounter
{


    //
    //Checks output data is roughly integer/actually percent
    //
    public void checkInteger(ReadCounter counter)
    {
        float[][] dat=counter.GetResults();
        String[] colNams=counter.GetColnames();
        int numCols=counter.numCol;
        int numCells=counter.numCell;
        boolean allGood=true;
        
        colLoop: 
        for(int i=1;i<(numCols-1);i++)
        {
            String curColumn=colNams[i+1];
            boolean isInt=true;//true if should be integer
            if(curColumn.length()>7){
                isInt=!curColumn.substring(0,7).equals("percent");
            }

            for(int j=0;j<numCells;j++)
            {
                if(isInt)
                {
                    if(Math.abs(dat[j][i]-Math.round(dat[j][i])) > .01)
                    {
                        allGood=false;
                        print("Results not close to integer");
                        print(colNams[i]);
                        print(Float.toString(dat[j][i]));
                        continue colLoop;
                    }
                }
                else{
                    if(dat[j][i]>100 | dat[j][i]<0)
                    {
                        allGood=false;
                        print("Percent not between 0 and 100, test failed");
                        print(Float.toString(dat[j][i]));
                        continue colLoop;
                    }
                }
            }
        }
        if(allGood)
        {
            print("Percents check out, as do counts!"); 
        }

    }

    //
    //compares output to metric output from CellRanger. Note do not expect perfect alignment, but most should be close.
    //
    public void compareToMetricCSV(ReadCounter counter,String metricFile)
    {
        String[] colsUseMetric={"Number of Reads","Reads Mapped Confidently to Intergenic Regions",
        "Reads Mapped Confidently to Intronic Regions","Reads Mapped Confidently to Exonic Regions",
        "Reads Mapped Antisense to Gene"};//columns to use in metric file
        
        String[] colsUseSC={"total","intergenic",
        "intronic","exonic","antisense"};//columns to use in counter object
        int numColsPrint=5;
        float[] vals_csv=loadMetricCSV(metricFile,colsUseMetric,numColsPrint);
        float[] vals_counter=aggregateCounterVals(counter,colsUseSC,numColsPrint);
        for(int i=0;i<numColsPrint;i++)
        {
            print("From metric csv, "+colsUseMetric[i]+": "+ Float.toString(vals_csv[i]));
            print("From counter, "+colsUseSC[i]+": "+ Float.toString(vals_counter[i]));
            print(" ");
        }
        
    }


    //extracts QC information from metric csv from CellRanger
    public float[] loadMetricCSV(String metricFile,String[] colsUseMetric,int numCol)
    {
        float[] ret=new float[numCol];
        Scanner s=null;
        try{
            s = new Scanner(new File(metricFile));
        }catch(Exception e){
            print("Exception with reading in metric file!");
        }
        String line=s.nextLine();
        String[] cols=line.split(","); //Metric files column names
        int numMetCols=1+(int)line.chars().filter(ch -> ch == ',').count();//number of columns in CSV


        //String line=s.nextLine().replace("%,","_").replace("\",","_");
        line=s.nextLine().replace("%","");
        //String[] vals=line.split("_"); //metric file values
        String[] vals=parseCSVLine(line,numMetCols);
        
        
        for(int i=0;i<numCol;i++)
        {
            int k=0;
            for(k=0;k<numMetCols;k++)
            {
                if(colsUseMetric[i].equals(cols[k].replace("\"","")))
                {
                    break;
                }
                
            }
            if(k==numMetCols)
            {
                print("Could not find column "+colsUseMetric[i]);
                continue;
            }
           
            ret[i]=Float.parseFloat(vals[k].replace("\"","").replace(",","").replace("%",""));

        }

        return(ret);
    }

    //splits apart lines of the metrics CSV file
    public String[] parseCSVLine(String line,int numMetCols)
    {
        boolean inQuotes=false; //if currently in quotations or not
        int lenString=line.length(); //length of string
        char[] charList=line.toCharArray(); //string to char array for easy modification

        for(int i=0;i<lenString;i++)
        {
            char curChar=charList[i];
            if(curChar=='\"')
            {
                if(inQuotes)
                {
                    inQuotes=false;
                }
                else{
                    inQuotes=true;
                }
                continue;
            }
            if(inQuotes)
            {
                continue;
            }
            if(curChar==',')
            {
                charList[i]='_';
            }

        }
        String lineNew=new String(charList);
        String[] ret=lineNew.split("_");
        return(ret);
    }

    //Aggregates the output of the counter class by column
    public float[] aggregateCounterVals(ReadCounter counter,String[] colsUseMetric,int numCol)
    {
        float[] ret=new float[numCol];

        int numCells=counter.numCell;

        float total=0; //Total number of reads
        for(int j=0;j<numCells;j++)
        {
            total=total+counter.CellQC[j][(counter.numCol-2)];
        }

        for(int i=0;i<numCol;i++)
        {
            int k=1;
            while(k<counter.numCol)
            {
                if(colsUseMetric[i]==counter.colNames[k])
                {
                    break;
                }
                k=k+1;
            }
            if(k==counter.numCol)
            {   
                print("Failed to find "+colsUseMetric[i]);
                ret[i]=-1;
                continue;
            }

            for(int j=0;j<numCells;j++)
            {
                ret[i]=ret[i]+counter.CellQC[j][(k-1)];
            }
            if(colsUseMetric[i]!="total" & colsUseMetric[i]!="percent_qual_cbc" & colsUseMetric[i]!="percent_qual_umi")
            {
                ret[i]=100*ret[i]/total;
            }
            
            
        }


        return ret;
    }

    public static void print(String output)
    {
        System.out.println(output);
    }

}
