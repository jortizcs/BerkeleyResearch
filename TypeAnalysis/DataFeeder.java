import java.nio.*;
import java.nio.channels.*;
import java.io.*;
import java.util.*;
import java.util.concurrent.*;
import java.text.*;

//hmm library
import be.ac.ulg.montefiore.run.jahmm.*;
import be.ac.ulg.montefiore.run.jahmm.toolbox.*;
import be.ac.ulg.montefiore.run.jahmm.learn.*;
import be.ac.ulg.montefiore.run.jahmm.draw.*;

import java.lang.*;



public class DataFeeder{
    private static File f;
    //file info
    private static final String filePath = null;
    private static RandomAccessFile file = null;
    private static MappedByteBuffer memoryBuffer = null;
    private static int count = 1024*1024*10; //10 MB
    private static int pos = 0;
    private static int readings = 100;

    //  [streamName, Hmm]
    private static ConcurrentHashMap<String, Hmm<ObservationReal>> modelMap=null;

    public static void main(String[] args){
        NativeClusterer nc = new NativeClusterer();
        try {
            DataFeeder feeder = new DataFeeder("datum");
            Matrix tempMatrix = feeder.analyze("DB1.txt", 1400);
            Matrix tempMatrix2 = feeder.analyze("DB1.txt",1400);
            Matrix tempMatrix3 = feeder.analyze("DB1.txt",1400);
            nc.addObservation("temperature", tempMatrix);
            nc.addObservation("temperature", tempMatrix2);
            nc.addObservation("temperature", tempMatrix3);
            feeder.reboot(null);
            Matrix rhMatrix = feeder.analyze("RH1.txt", 1400);
            Matrix rhMatrix2 = feeder.analyze("RH1.txt", 1400);
            Matrix rhMatrix3 = feeder.analyze("RH1.txt", 1400);
            nc.addObservation("rh", rhMatrix);
            nc.addObservation("rh", rhMatrix2);
            nc.addObservation("rh", rhMatrix3);

            //test classifier
            feeder.reboot(null);
            Matrix testMatrix = feeder.analyze("RH1.txt", 1400);
            Matrix testMatrix2 = feeder.analyze("RH1.txt", 1400);
            Matrix testMatrix3 = feeder.analyze("RH1.txt", 1400);
            ArrayList<Matrix> dps = new ArrayList<Matrix>();
            dps.add(testMatrix); dps.add(testMatrix2);dps.add(testMatrix3);

            feeder.reboot(null);
            Matrix testMatrix4 = feeder.analyze("DB1.txt", 1400);
            Matrix testMatrix5 = feeder.analyze("DB1.txt", 1400);
            Matrix testMatrix6 = feeder.analyze("DB1.txt", 1400);
            ArrayList<Matrix> dps2 = new ArrayList<Matrix>();
            dps2.add(testMatrix4); dps2.add(testMatrix5);dps2.add(testMatrix6);


            System.out.println("Real Class: rh, Classification: " + nc.classify(dps));
            System.out.println("Real Class: temperature, Classification: " + nc.classify(dps2));
        } catch(Exception e){
            e.printStackTrace();
        }
    }

    public DataFeeder(String dataFilesDirPath) throws Exception{
        f = new File(dataFilesDirPath);
        if(!f.isDirectory())
            throw new Exception (dataFilesDirPath + " must be a directory");
        modelMap = new ConcurrentHashMap<String, Hmm<ObservationReal>>();
    }

    public void reboot(String dataFilesDirPath) throws Exception{
        if(dataFilesDirPath!=null)
            f = new File(dataFilesDirPath);
        if(!f.isDirectory())
            throw new Exception (dataFilesDirPath + " must be a directory");
        modelMap = new ConcurrentHashMap<String, Hmm<ObservationReal>>();
    }

    public Matrix analyze(String streamName, int sampleSize){
        try {
            String filePath = f.getCanonicalPath() + "/" + streamName;
            file = new RandomAccessFile(filePath, "r");
            memoryBuffer = file.getChannel().map(FileChannel.MapMode.READ_ONLY, 0, count);

            System.out.println("Populating sample vector");
            ArrayList<ObservationReal> obs = new ArrayList<ObservationReal>();
            for(int i=0; i<sampleSize; i++)
                obs.add(new ObservationReal(nextVal()));
            ArrayList<ArrayList<ObservationReal>> samples = 
                new ArrayList<ArrayList<ObservationReal>>();
            samples.add(obs);

            //learn and record the model for 1-10 states
            int maxStates = 8;
            ArrayList<Hmm<ObservationReal>> models = new ArrayList<Hmm<ObservationReal>>();
            for(int i=1; i<=maxStates; i++){
                System.out.print("Learning Models:"); 
                System.out.print(i + " of " + maxStates + " states, samples size=" + samples.get(0).size());
                KMeansLearner<ObservationReal> kml = 
                    new KMeansLearner<ObservationReal>(1, new OpdfGaussianFactory(), samples);
                try {
                    System.out.println("learning for " + i + " states");
                    Hmm<ObservationReal> model = kml.learn();
                    if(i==2 || i==4)
                        (new GenericHmmDrawerDot()).write(model, 
                                streamName + "_learntHmm_"+ i +".dot");
                    System.out.println("done");
                    models.add(i-1,model);
                } catch(Exception e){
                    e.printStackTrace();
                }
                System.out.print("\r");
            }

            //choose model with least error
            System.out.println("\nPicking best model");
            int k = 8*1000;
            ArrayList<ObservationReal> nextK = new ArrayList<ObservationReal>();
            Hmm<ObservationReal> bestModel = null;
            double minAvgError = Double.MAX_VALUE;
            for(int j=0; j<k; j++)
                nextK.add(new ObservationReal(nextVal()));
            for(int i=0; i<maxStates; i++){
                Hmm<ObservationReal> thisModel = models.get(i);
                if(thisModel!=null){
                    MarkovGenerator gen = new MarkovGenerator(thisModel);
                    ArrayList<ObservationReal> sequence = 
                        new ArrayList<ObservationReal>(gen.observationSequence(k));
                    double avgError = 0L;
                    double totError = 0L;
                    for(int j=0; j<k; j++){
                        double abserr = sequence.get(j).value - nextK.get(j).value;
                        totError += Math.pow(abserr,2);
                    }
                    avgError = totError/k;
                    System.out.println("states=" + (i+1) + ", error=" + avgError);
                    if(avgError<minAvgError){
                        minAvgError = avgError;
                        bestModel = thisModel;
                    }
                }
            }

            //record the best model
            modelMap.put(streamName, bestModel);

            //close the stream data file
            file.close();

            //done
            System.out.println("Done");
            System.out.println(bestModel);
            Matrix m = new Matrix(9,9);
            System.out.print("[\n");
            for(int i=0; i<9; i++){
                for(int j=0; j<9; j++){
                    DecimalFormat f = new DecimalFormat("#0.000");
                    if(i<bestModel.nbStates()-1 && j<bestModel.nbStates()-1){
                        System.out.print(" " + f.format(bestModel.getAij(i, j)));
                        m.put(i,j,bestModel.getAij(i, j));
                    }else{
                        System.out.print(" 0.000");
                        m.put(i,j,0);
                    }
                }
                System.out.println();
            }
            System.out.println("]");
            return m;

        } catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
        return null;
    }

    public double nextVal(){
        StringBuffer line = new StringBuffer();
		int cnt = 0;
		try {
			char c = (char)memoryBuffer.get();
			while(c!='\0' && c!='\n'){
				line.append(c);
				c = (char)memoryBuffer.get();
			}
			//System.out.println(line);
			StringTokenizer tokenizer = new StringTokenizer(line.toString(), "\t");
			tokenizer.nextToken();
			String v = tokenizer.nextToken();
			if(v!=null)
				return Double.parseDouble(v);
		} catch(Exception e){
			e.printStackTrace();
		}
		return 0.0;
    }

}
