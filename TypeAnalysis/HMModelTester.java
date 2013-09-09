import java.io.*;
import java.nio.*;
import java.nio.channels.*;
import java.util.*;
import java.util.concurrent.locks.*;

//hmm library
import be.ac.ulg.montefiore.run.jahmm.*;
import be.ac.ulg.montefiore.run.jahmm.toolbox.*;
import be.ac.ulg.montefiore.run.jahmm.learn.*;

import java.lang.*;

public class HMModelTester{
    private int nsamples = 6000;
    private ArrayList<ArrayList<ObservationReal>> samples;
    private int index = 0;
    private int totalReadings = 0;

    private long lastReceivedTs = 0L;
    private long totalPeriods = 0L;
    private long avgperiod = 0L;
    private long lastHmmRun = 0L;

    private static ByteBuffer rbuffer = null;

    //hmm-related
    private int nbstates = 3;
    private static Hmm<ObservationReal> hmm = null;
    private static Hmm<ObservationReal> hmm2 = null;
    private MarkovGenerator gen = null;

    //file info
    private static final String filePath = "DB1.txt";
    private static RandomAccessFile file = null;
    private static MappedByteBuffer memoryBuffer = null;
    private static int count = 1024*1024*10; //10 MB
    private static int pos = 0;
    private static int readings = 100;

    public HMModelTester(){
        samples = new ArrayList<ArrayList<ObservationReal>>(readings);
    }

    public static void main(String[] args){
        HMModelTester hmodel = new HMModelTester();
        rbuffer =  ByteBuffer.allocate(16*readings);
        long now  = System.currentTimeMillis()-(readings * 1000*60*15);

        //open data file and read it
        try {
            file = new RandomAccessFile(filePath, "r");
            memoryBuffer = file.getChannel().map(FileChannel.MapMode.READ_ONLY, 0, count);
            long ts = now;
            for(int j=0; j<2; j++){
                int i=0;
                while(i<readings){
                    rbuffer.putLong(ts);
                    rbuffer.putDouble(nextTempVal());
                    ts+=(1000*60*15);
                    hmodel.dataReceived();
                    i+=1;
                }
                rbuffer.clear();
                if(hmm2==null)
                    hmm2 = hmm.clone();
            }
        } catch(Exception e){
            e.printStackTrace();
        } finally {
            try{
                file.close();
            } catch(Exception e){
            }
            KullbackLeiblerDistanceCalculator klc = new KullbackLeiblerDistanceCalculator();
            System.out.println("Distance:" + klc.distance(hmm,hmm2));
        }
    }

	public static double nextTempVal(){
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
    

    public void setNoStates(int nstates){
        nbstates = nstates;
    }

    public int getNoStates(){
        return nbstates;
    }

    private void reset(){
        hmm = null;
        index = 0;
        samples.clear();
    }

    public boolean ready(){
        return hmm!=null;
    }

    public double valueAt(long ts){
       if(lastHmmRun>0){
           if(avgperiod>=0){
               int seqLength =(int) ((double)(ts-lastHmmRun)/avgperiod);
               if(seqLength>1000)
                   seqLength=1;
               return ((ObservationReal)gen.observationSequence(seqLength).get(seqLength-1)).value;
           } 
       }
       if(samples.size()>0)
           return samples.get(samples.size()-1).get(0).value;
       return 0.0;
    }

    private static ArrayList<ObservationReal> obs = new ArrayList<ObservationReal>(readings);
    public void dataReceived(){
        try{
            int max = readings -2;
            System.out.print(index + "/" + max); Thread.sleep(10);
            System.out.print('\r');
        } catch(Exception e){}
        if(index>=readings-2){
            samples.add(0, obs);
            refit();
            index = 0;
            samples.clear();
        }

        //ArrayList<ObservationReal> obs = new ArrayList<ObservationReal>(1);
        //obs.add(new ObservationReal(rbuffer.getLong(rbuffer.position()-16)));
        obs.add(new ObservationReal(rbuffer.getDouble(rbuffer.position()-8)));
        //samples.add(index, obs);
        index+=1;
    }

    public void refit(){
        long before = System.currentTimeMillis();
        for(int i=2; i<=nbstates; i++){
            try {
                KMeansLearner<ObservationReal> kml = 
                    new KMeansLearner<ObservationReal>(i, 
                            new OpdfGaussianFactory(), samples);
                if(i==2)
                    hmm = kml.learn();
                else
                    kml.learn();
            } catch(Exception e){
                System.out.println("nbstates=" + i + " fails");
                //e.printStackTrace();
            }
        }
        gen = new MarkovGenerator(hmm);
        long time_ = System.currentTimeMillis()-before;
        System.out.println("HMM Computed::time=" + time_);
        System.out.println(hmm);
        //System.exit(0);
    }
}
