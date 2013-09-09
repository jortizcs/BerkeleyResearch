import java.util.*;
import java.text.*;
import java.io.*;

public class DataProcessor{

    //public static final String tsFormat = "yyyy-MM-DD'T'HH:mm:ss'-07:00'";
    public static final String tsFormat = "yyyy-MM-DD HH:mm:ss'+09'";

    public DataProcessor(){}

    public static void main(String[] args){
        try {
            DataProcessor processor = new DataProcessor();
            File[] dataFiles = processor.fileList("../data/week");
            for(int i=0; i<dataFiles.length; i++){
                if(!dataFiles[i].isDirectory()){
                    System.out.println(dataFiles[i]);
                    StringTokenizer tokenizer = new StringTokenizer(dataFiles[i].getCanonicalPath(),"/");
                    int numTokens = tokenizer.countTokens();
                    int cnt = 0;
                    StringBuffer dirPath = new StringBuffer();
                    while(cnt<numTokens-1){
                        cnt+=1;
                        dirPath.append(tokenizer.nextToken()).append("/");
                    }
                    String name = tokenizer.nextToken();
                    tokenizer = new StringTokenizer(name, ".");
                    name = tokenizer.nextToken();
                    name  = new StringBuffer().append("/").append(dirPath.toString()).
                       append("processed/").append(name).append(".processed.csv").
                       toString();
                    System.out.println(name);
                    processor.parseFile(dataFiles[i], name, "\t", DataProcessor.tsFormat);
                }
            }
        } catch(Exception e){
            e.printStackTrace();
        }
    }

    public File[] fileList(String rootDirPath){
        ArrayList<File> flist = new ArrayList<File>();
        File rootDirFile = new File(rootDirPath);
        if(rootDirFile.isDirectory())
            return rootDirFile.listFiles();
        return null;
    }

    public void parseFile(File dataFile, String destPath, String delim, String timeStampFormat){
        try {
            BufferedReader fileReader = new BufferedReader(new FileReader(dataFile));
            String line = null;
            SimpleDateFormat formatter = new SimpleDateFormat(timeStampFormat);
            BufferedWriter writer = new BufferedWriter(new FileWriter(destPath, true));
            while((line=fileReader.readLine())!=null){
                StringTokenizer tokenizer = new StringTokenizer(line, delim);
                if(tokenizer.countTokens()>=2){
                    try {
                        String ts = tokenizer.nextToken();
                        String v = tokenizer.nextToken();
                        //System.out.println("t="+ts + ", v=" + v);
                        String value = new Double(Double.parseDouble(v)).toString();

                        Date tsDate = null;
                        if(ts!=null && (tsDate = formatter.parse(ts, new ParsePosition(0)))!=null){
                            StringBuffer lineBuf = new StringBuffer().
                                append(new Long(tsDate.getTime()).toString()).
                                append(",").append(value).append("\n");
                            String lineStr = lineBuf.toString();
                            writer.write(lineStr, 0, lineStr.length());
                            //System.out.print(lineStr);
                        }
                    } catch(Exception e){
                        e.printStackTrace();
                    }
                }
            }
            writer.flush();
        } catch(Exception e){
            e.printStackTrace();
        }
    }
}
