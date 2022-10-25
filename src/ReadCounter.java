package SingleCellQC;
import java.io.FileWriter;
import java.io.*;
import htsjdk.samtools.*;
import java.util.*;
import java.util.Map.Entry;
import java.util.Map;
import java.lang.*;
import java.time.Duration;
import java.time.Instant;
import java.util.zip.GZIPInputStream;

//////////////////////////////////////////////////////////
//// This class functions to read in and store info about the number of reads per cell
//// in certain classes (intronic, intergenic, etc) in a CellRanger run. Built to run
//// with CellRanger v6, should work with other versions as well (though some values
//// might be off). Note assumes not using a pre-mrna reference (if you are exonic and 
//// intronic read numbers will be off) though works with the include introns flag.
//////////////////////////////////////////////////////////
public class ReadCounter
{

    //Column names
    protected String[] colNames={"CBC","antisense","intergenic","intronic","exonic","multi","unmapped","highConf","polyA","TSO","total"}; //Column names in output

    //location columns
    protected final int col_anti=0; //column with count of antisense reads
    protected final int col_intergenic=1; //column with count of intergenic reads
    protected final int col_intronic=2; //column with count of intronic reads
    protected final int col_exonic=3; //column with count of exonic reads
    protected final int col_multi=4; //column with count of multimapped reads
    protected final int col_unmap=5; //column with count of unmapped reads
    protected final int col_hiconf=6; //column with count of reads with high confidence mappings to transcriptome
    protected final int col_polyA=7; //column with count of reads with polyA trimmed off (only counts aligned reads)
    protected final int col_TSO=8; //column with count of reads with TSO trimmed off (only counts aligned reads)
    protected final int col_tot=9; //column with count of all reads
    protected final int numCol=11; //Number of columns, including CBC (so the value in col_tot plus 2 if col_tot is last column)

    //Files
    protected File bamFile; //The CellRanger bam file
    protected File cellFile; //A text file with the cell barcodes being read in, assume not gzipped (will implement later)
    protected File outfile; //The name of the file to write to, if exists will be overwritten

    //Data structures
    protected ArrayList<String> cells; //list of CBCs
    protected int numCell; //Number of CBC
    protected HashMap<String, Integer> Cell2Pos; //Maps from cell barcode to position in cells
    protected float[][] CellQC; //A 2 dimensional array, rows are cells, columns correspond to the different quantities of interest

    
    ////////////////////
    ////Initialize object for counting
    ////bamFile: The name of the bam file from CellRanger
    ////cellFile: The name of the file containing CBCs to look for
    ////outfile: The name of the file used for saving the results
    ////gzip: true if cell name file is gzipped
    /////////////////////
    public ReadCounter(String bamFile,String cellFile,String outfile,boolean gzip)
    {
        this.bamFile=new File(bamFile); //The bam file to process
        this.cellFile=new File(cellFile); //the file with CBC info
        this.outfile=new File(outfile);//the output file
        print("Read in Cell Data");
        this.cells=new ArrayList<String>();
        try{
            Scanner s;
            if(!gzip)
            {
                s = new Scanner(this.cellFile);
            }
            else
            {
                s = new Scanner(new GZIPInputStream(new FileInputStream(this.cellFile)));
            }
            while(s.hasNext())
            {
                this.cells.add(s.next());
            }
           
            s.close();
        }catch(Exception e){
            print("Exception with reading in cell file!");
        }
        this.numCell=this.cells.size();

        this.Cell2Pos=new HashMap<String, Integer>(); //Maps CBC to positive in array
        for(int i=0;i<this.numCell;i++)
        {
            String cell=this.cells.get(i);
            
            this.Cell2Pos.put(cell,i);
        }

        

    }


    //Reads each alignment in the bam one by one and gets QC info
    public void ReadBam(boolean verbose)
    {
        print("Read in data!");
        SamReader sr = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(this.bamFile);
        SAMRecordIterator r = sr.iterator();


        this.CellQC=new float[this.numCell][numCol-1]; //Stores the QC information we care about

        int readNum=0;
        int readCounted=0;


        
        Instant inst1 = Instant.now();
        while(r.hasNext()) {

            SAMRecord read=r.next(); //Current Read
            readNum=readNum+1;
            if(readNum % 1000000==0 & verbose)
            {
                print(String.valueOf(readNum));
                print("");
            }
            
            //To speed up testing
            //if(readNum>100000)
            //{
            //    break;
            //}
            ///////
        
            String cbc=""; //current cell barcode
            try{
                cbc=read.getStringAttribute("CB");
            }catch(Exception e){
                cbc="NotCell";
            }

            

            if(!this.Cell2Pos.containsKey(cbc)) 
            {
                continue;
            }

            readCounted=readCounted+1; //counts number of reads with CBC in list
            int pos=this.Cell2Pos.get(cbc); //row this cbc appears in


            int numMapping; //Number position read maps to
            try{
                numMapping=read.getIntegerAttribute("NH");
            }
            catch(Exception e)
            {
                numMapping=0;
            }


            if(numMapping<1)
            {
                this.CellQC[pos][this.col_unmap]=this.CellQC[pos][this.col_unmap]+1;
                this.CellQC[pos][this.col_tot]=this.CellQC[pos][this.col_tot]+1;
                continue;
            }

            int numPolyA; //Number of bases trimmed for polyA
            int numTSO; //Number of bases timmed for TSO
            try{
                numPolyA=read.getIntegerAttribute("pa");
            }
            catch(Exception e)
            {
                numPolyA=0;
            }

            float numMapping_float=numMapping;

            

            if(numPolyA>0)
            {
                this.CellQC[pos][this.col_polyA]=this.CellQC[pos][this.col_polyA]+1/numMapping_float;
            }

            try{
                numTSO=read.getIntegerAttribute("ts");
            }
            catch(Exception e)
            {
                numTSO=0;
            }
            
            if(numTSO>0)
            {
                this.CellQC[pos][this.col_TSO]=this.CellQC[pos][this.col_TSO]+1/numMapping_float;
            }
            


            if(numMapping>1)
            {
                this.CellQC[pos][this.col_tot]=this.CellQC[pos][this.col_tot]+1/numMapping_float;
                this.CellQC[pos][this.col_multi]=this.CellQC[pos][this.col_multi]+1/numMapping_float;
                continue;
            }

            this.CellQC[pos][this.col_tot]=this.CellQC[pos][this.col_tot]+1;

            char readType; //If intergenic, intornic, or exonic
            try{
                readType=read.getCharacterAttribute("RE");
            }
            catch(Exception e)
            {
                continue;
            }

            if(readType=='E')
            {
                this.CellQC[pos][this.col_exonic]=this.CellQC[pos][this.col_exonic]+1; //If exonic
            }

            if(readType=='N') 
            {
                this.CellQC[pos][this.col_intronic]=this.CellQC[pos][this.col_intronic]+1; //If intronic
            }

            if(readType=='I')
            {
                this.CellQC[pos][this.col_intergenic]=this.CellQC[pos][this.col_intergenic]+1; //If intergenic
            }

            int xf=0; //to check if confidentially mapped to transcriptome

            try{
                xf=read.getIntegerAttribute("xf");
            }
            catch(Exception e)
            {
                continue;
            }

            if(xf % 2==1)
            {
                this.CellQC[pos][this.col_hiconf]=this.CellQC[pos][this.col_hiconf]+1; //Adds if hi confidence map to transcriptome
            }

            
            String antisense;
            try
            {
                antisense=read.getStringAttribute("AN");
            }
            catch(Exception e)
            {
                continue;
            }

            String gene;

            try
            {
                gene=read.getStringAttribute("TX");
            }
            catch(Exception e)
            {
                continue;
            }

            if(gene==null & !(antisense==null))
            {
                this.CellQC[pos][this.col_anti]=this.CellQC[pos][this.col_anti]+1; //Adds for antisense
            }


        }
        
        Instant inst2 = Instant.now(); 
        
        print(Duration.between(inst1, inst2).toString());
        
        print("Total number of alignments: "+String.valueOf(readNum));

    }

    public void SaveQC()
    {

        print("Save file");


        try
        {
            if(!this.outfile.exists())
            {
                this.outfile.createNewFile();
            }
            FileWriter fw = new FileWriter(this.outfile.getAbsoluteFile());
            BufferedWriter bw = new BufferedWriter(fw);
            for(int j=0;j<numCol;j++)
            {
                String colnam=colNames[j];
                if(j>0)
                {
                    colnam="\t"+colnam;
                }
                bw.write(colnam);
            }
            bw.newLine();

            for(int i=0;i<this.numCell;i++)
            {
                String cell=this.cells.get(i);
                float[] quant=this.CellQC[i];
                bw.write(cell,0,cell.length()); //write Cell name
                for(int j=0;j<numCol-1;j++)
                {
                    bw.write("\t"+String.valueOf(Math.round(quant[j])));
                }
                bw.newLine();
            }
            bw.close();
        }
        catch(Exception e)
        {
            print("Exception writing file!");
        }


    }


    public static void print(String output)
    {
        System.out.println(output);
    }

}
