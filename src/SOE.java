import java.util.Arrays;

public class SOE {
    GlobalData globalData = new GlobalData();
    public double [][] HG = new double[globalData.getnN()][globalData.getnN()];
    public double [][] CG = new double[globalData.getnN()][globalData.getnN()];
    public double [] PG = new double[globalData.getnN()];
    public double [] t1 = new double[globalData.getnN()];

    public SOE(){};
}
