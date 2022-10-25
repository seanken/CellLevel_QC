package SingleCellQC;
import java.util.*;
import java.lang.*;


public class CellLevel_QC
{
    public static void main(String[] args)
    {
    ReadCounter counter=new ReadCounter(args[0],args[1],args[2]);
    counter.ReadBam();
    counter.SaveQC();

    }


    public static void print(String output)
    {
        System.out.println(output);
    }

}
