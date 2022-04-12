import java.util.Arrays;

public class Element {

    public int[] id;
    public double [][] Hl = new double[4][4];
    public double [][] Cl = new double[4][4];
    public double [] Pl = new double[4];

    public Element(int[] id) {
        this.id = id;
    }
    public Element() {};

    @Override
    public String toString() {
        return "Element{" +
                "id=" + Arrays.toString(id) +
                '}';
    }
}
