import java.io.File;
import java.io.FileNotFoundException;
import java.util.Scanner;

public class GlobalData {

    private double H;               //wysokość
    private double W;               //szerokość
    private int nH;                 //liczba węzłów po wysokości
    private int nW;                 //liczba węzłów po szerokości
    private int nN;                 //liczba węzłów
    private int nE;                 //liczba elementów
    private int schematc;           //schematcalkowania
    private double k;               //przewodność cieplna
    private double c;               //ciepło właściwe
    private double ro;              //gęstość
    private double alfa;            //alfa
    private double totoczenia;      //temperatura otoczenia
    private double deltat;          //krok czasowy symulacji
    private double t0;              //temperatura początkowa
    public double tsymulacji;       //czas symulacji

    public GlobalData() {
        try {
            File file = new File("dane.txt");
            Scanner scanner = new Scanner(file);
            this.H = scanner.nextDouble();
            this.W = scanner.nextDouble();
            this.nH = scanner.nextInt();
            this.nW = scanner.nextInt();
            this.schematc = scanner.nextInt();
            this.k=scanner.nextDouble();
            this.c=scanner.nextDouble();
            this.ro=scanner.nextDouble();
            this.alfa=scanner.nextDouble();
            this.totoczenia=scanner.nextDouble();
            this.deltat=scanner.nextDouble();
            this.t0=scanner.nextDouble();
            this.tsymulacji=scanner.nextDouble();

            scanner.close();
        } catch (FileNotFoundException f) {
            System.out.println("Blad pliku!");
            f.printStackTrace();
        }
        this.nN = nH * nW;
        this.nE = (nH - 1) * (nW - 1);
    }
    public double getH() {
        return H;
    }
    public double getW() {
        return W;
    }
    public int getnH() {
        return nH;
    }
    public int getnW() {
        return nW;
    }
    public int getnN() {
        return nN;
    }
    public int getnE() { return nE; }
    public int getSchematc() { return schematc; }
    public double getK() { return k; }
    public double getC() { return c; }
    public double getRo() { return ro; }
    public double getAlfa() { return alfa; }
    public double getTotoczenia() { return totoczenia; }
    public double getDeltat() { return deltat; }
    public double getT0() { return t0; }
    public double getTsymulacji() { return tsymulacji; }
}
