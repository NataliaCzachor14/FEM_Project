public class Node {
    public double x;
    public double y;
    public double BC;
    public double t0;

    public Node(double x, double y, double BC, double t0) {
        this.x = x;
        this.y = y;
        this.BC = BC;
        this.t0 = t0;
    }

    public Node() {};

    public double getX() {
        return x;
    }

    public double getY() {
        return y;
    }

    public double getBC() { return BC; }

    public double getT0() { return t0; }

    public void setT0(double t0) { this.t0 = t0; }

    @Override
    public String toString() {
        return "Node{" +
                "x=" + x +
                ", y=" + y +
                '}';
    }
}
