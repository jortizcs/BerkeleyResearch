import java.util.*;
import java.util.concurrent.*;

public class Matrix{
    private ConcurrentHashMap<String, Double> internal;
    public int rows = 0;
    public int cols = 0;

    public Matrix(int nrows, int ncols){
        setDims(nrows, ncols);
        internal = new ConcurrentHashMap<String, Double>();
    }

    public static void main(String[] args){
        Matrix testMatrix = new Matrix(9,9);
        testMatrix.put(0,0,1);
        System.out.println("testMatrix(0,0)="+testMatrix.get(0,0));
        System.out.println(testMatrix.internal);
    }

    public void setDims(int nrows, int ncols){
        if(nrows>=0 && ncols >=0){
            rows = nrows;
            cols = ncols;
        }
    }

    public void put(int i, int j, double val){
        if(i>=0 && j>=0 && i<rows && j<cols){
            String key = new StringBuffer().append(i).append("+").append(j).toString();
            //System.out.println("key_put(key=" + key + ")");
            if(!internal.containsKey(key))
                internal.put(key, new Double(val));
            else
                internal.replace(key, new Double(val));
        }
    }
    
    public double get(int i, int j){
        if(i>=0 && j>=0 && i<rows && j<cols){
            String key = new StringBuffer().append(i).append("+").append(j).toString();
            if(internal.containsKey(key))
                return internal.get(key).doubleValue();
            else
                return 0.0;
        }
        return Double.NaN;
    }

    public int nrows(){
        return rows;
    }

    public int ncols(){
        return cols;
    }
}
