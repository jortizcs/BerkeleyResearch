import java.util.*;
import java.text.*;
import java.io.*;

public class OpenBmsDataProcessor{

    public static final String tsFormat = "yyyy-MM-DD'T'HH:mm:ss'-07:00'";

    public OpenBmsDataProcessor(){}

    public static void main(String[] args){
        String dataDir = "../data";
        String procDir = null;
        String filePath = null;

        // get the arguments
        int input=0;
        while (input<args.length){
            if(input<args.length-1){
                if(args[input].equals("--data") )
                    dataDir = args[input+1];
                else if (args[input].equals("--processed"))
                    procDir = args[input+1];
                else if (args[input].equals("--file"))
                    filePath = args[input+1];
            }
            input+=1;
        }

        if(filePath!=null){
            try {
                OpenBmsDataProcessor processor = new OpenBmsDataProcessor();
                File dataFile = new File(filePath);
                if(!dataFile.isDirectory()){
                    System.out.println(dataFile);
                    StringTokenizer tokenizer = new StringTokenizer(dataFile.getCanonicalPath(),"/");
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
                       append("processed/").append(name).append(".processed.csv").toString();
                    if(procDir!=null)
                        name  = new StringBuffer().append(procDir).append(name).append(".processed.csv")                                    .toString();
                    System.out.println(name);
                    processor.parseFile(dataFile, name, ",", OpenBmsDataProcessor.tsFormat);
                }
            } catch(Exception e){
                e.printStackTrace();
            }
                
        } else {
            // parse it, process it
            try {
                OpenBmsDataProcessor processor = new OpenBmsDataProcessor();
                File[] dataFiles = processor.fileList(dataDir);
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
                           append("processed/").append(name).append(".processed.csv").toString();
                        if(procDir!=null)
                            name  = new StringBuffer().append(procDir).append(name).append(".processed.csv")                                    .toString();
                        System.out.println(name);
                        processor.parseFile(dataFiles[i], name, ",", OpenBmsDataProcessor.tsFormat);
                    }
                }
            } catch(Exception e){
                e.printStackTrace();
            }
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
                    String ts = tokenizer.nextToken();
                    String value = new Double(Double.parseDouble(tokenizer.nextToken())).toString();

                    Date tsDate = null;
                    if(ts!=null && (tsDate = formatter.parse(ts, new ParsePosition(0)))!=null){
                        StringBuffer lineBuf = new StringBuffer().
                            append(new Long(tsDate.getTime()).toString()).
                            append(",").append(value).append("\n");
                        String lineStr = lineBuf.toString();
                        writer.write(lineStr, 0, lineStr.length());
                        //System.out.print(lineStr);
                    }
                }
            }
            writer.flush();
        } catch(Exception e){
            e.printStackTrace();
        }
    }
}
