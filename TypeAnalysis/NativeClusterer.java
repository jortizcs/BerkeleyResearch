import java.util.*;
import java.util.concurrent.*;

public class NativeClusterer{
    private ConcurrentHashMap<String, Matrix> centroids=null;
    private ConcurrentHashMap<String, Integer> pointsCnt=null;
   
    public NativeClusterer(){
        centroids = new ConcurrentHashMap<String, Matrix>();
        pointsCnt = new ConcurrentHashMap<String, Integer>();
    }

    public void addObservation(String name, Matrix m){
        if(centroids.containsKey(name)){
            Matrix oldCentroid = centroids.get(name);
            Matrix newCentroid = new Matrix(oldCentroid.nrows(), oldCentroid.ncols());
            int totalPoints = pointsCnt.get(name);
            for(int i=0; i<m.nrows(); i++){
                for(int j=0; j<m.ncols(); j++){
                    double v =  m.get(i,j) + 
                            (oldCentroid.get(i,j)*totalPoints)/(totalPoints+1);
                    newCentroid.put(i, j, v);
                }
            }
            totalPoints +=1;
            pointsCnt.put(name, new Integer(totalPoints));
            centroids.put(name, newCentroid);
        } else {
            pointsCnt.put(name, new Integer(1));
            centroids.put(name, m);
        }
    }

    public String classify(Matrix m){
        Iterator<String> keys = centroids.keySet().iterator();
        String mClass = null;
        double minDist = Double.MAX_VALUE;
        while(keys.hasNext()){
            String cls = keys.next();
            double dist = getDistance(m, centroids.get(cls));
            System.out.println("class:" + cls + ", dist:" + dist);
            if(dist<minDist){
                minDist = dist;
                mClass = cls;
            }
        }
        return mClass;
    }

    public String classify(List<Matrix> pts){
        //calculate the centroid
        Matrix centroid = null;
        for(int idx=0; idx<pts.size(); idx++){
            Matrix thisDatapt = pts.get(idx);
            if(centroid==null)
                centroid = new Matrix(thisDatapt.nrows(), thisDatapt.ncols());
            for(int i=0; i<thisDatapt.nrows(); i++){
                for(int j=0; j<thisDatapt.ncols(); j++){
                    double v =  thisDatapt.get(i,j) + 
                            (centroid.get(i,j)*i)/(i+1);
                    centroid.put(i, j, v);
                }
            }
        }

        //compare to the others
        Iterator<String> keys = centroids.keySet().iterator();
        String mClass = null;
        double minDist = Double.MAX_VALUE;
        while(keys.hasNext()){
            String cls = keys.next();
            double dist = getDistance(centroid, centroids.get(cls));
            System.out.println("class:" + cls + ", dist:" + dist);
            if(dist<minDist){
                minDist = dist;
                mClass = cls;
            }
        }
        return mClass;
    }

    private Matrix getCentroid(String name){
        return centroids.get(name);
    }

    private double getDistance(Matrix m1, Matrix m2){
        double dist = 0.0;
        System.out.println("nrows="+m1.nrows()+", ncols=" + m1.ncols());
        for(int i=0; i<m1.nrows(); i++){
            for(int j=0; j<m1.ncols(); j++){
                double diff = m1.get(i,j)-m2.get(i,j);
                dist += Math.pow(diff,2);
            }
        }
        return Math.sqrt(dist);
    }
}
