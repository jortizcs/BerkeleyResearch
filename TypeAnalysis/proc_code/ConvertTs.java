import java.util.*;
import java.text.*;
import java.io.*;

public class ConvertTs{
    public static void main(String[] args){
        String formatStr = null;
        String timeStr = null;
        for(int i=0; i<args.length; i++){
            //System.out.println("args[" + i + "]=" + args[i]);
            if(args[i].equals("--format") && i+1<args.length){
                formatStr = args[i+1];
            } else if(args[i].equals("--time") && i+1<args.length){
                timeStr = args[i+1];
            }
        }
        //System.out.println("formatStr=" + formatStr+ ", timeStr=" + timeStr);

        if(timeStr!=null && formatStr!=null){
            SimpleDateFormat formatter = new SimpleDateFormat(formatStr);
            Date tsDate = new Date(Long.parseLong(timeStr));
            StringBuffer buf = new StringBuffer();
            System.out.println(formatter.format(tsDate, buf, new FieldPosition(0)));
        }
    }
}
