import java.util.ArrayList;
import java.util.List;

public class Main {

    static GlobalData globalData = new GlobalData();
    static List<Node> nodeList = new ArrayList<>();
    static List<Element> elementList = new ArrayList<>();
    static SOE soe = new SOE();

    public static void main(String[] args) {

        //tworzenie nodeList
        double deltax = (globalData.getW()) / (globalData.getnW() - 1);
        double deltay = (globalData.getH()) / (globalData.getnH() - 1);

        for (int i = 0; i < globalData.getnW(); i++) {
            for (int j = 0; j < globalData.getnH(); j++) {
                if (j == 0 || j == (globalData.getnH() - 1) || i == 0 || i == (globalData.getnW() - 1))
                    nodeList.add(new Node(i * deltax, j * deltay, 1, globalData.getT0()));
                else {
                    nodeList.add(new Node(i * deltax, j * deltay, 0, globalData.getT0()));
                }
            }
        }

        //tworzenie elementList
        for (int i = 1; i < globalData.getnE() + (globalData.getnW() - 1); i++) {

            if (i % globalData.getnH() == 0) {
                continue;
            }
            Element el = new Element(new int[4]);
            el.id[0] = i;
            el.id[1] = el.id[0] + globalData.getnH();
            el.id[2] = el.id[1] + 1;
            el.id[3] = el.id[0] + 1;
            elementList.add(el);
        }

        //wypisywanie siatki
        /*for (int i = 0; i < globalData.getnE(); i++) {
            System.out.println("Element: " + (i + 1));
            System.out.println("Nodes: ");
            for (int j = 0; j < 4; j++) {
                System.out.print(elementList.get(i).id[j] + ": ");
                System.out.print("X: " + nodeList.get(elementList.get(i).id[j] - 1).getX() + ", ");
                System.out.print("Y: " + nodeList.get(elementList.get(i).id[j] - 1).getY()+ ", ");
                System.out.println("BC: " + nodeList.get(elementList.get(i).id[j] - 1).getBC()+ ", ");
            }
        }*/
        System.out.println();

        Elem4 elem4 = new Elem4(globalData.getSchematc());


        double[][] temp_H;     //H
        double[][] temp_C;     //C
        double[][] temp_HBC;   //HBC
        double[] temp_P;       //P

        for (int w = 0; w < (globalData.getTsymulacji() / globalData.getDeltat()); w++) {

            //System.out.println("\nIteracja: " + w);

            for (int i = 0; i < globalData.getnE(); i++) {
                for (int j = 0; j < 4; j++) {
                    for (int v = 0; v < 4; v++) {
                        elementList.get(i).Hl[j][v] = 0;
                        elementList.get(i).Cl[j][v] = 0;
                        elementList.get(i).Pl[j] = 0;
                    }
                }
            }

            for (int a = 0; a < globalData.getnN(); a++) {
                for (int b = 0; b < globalData.getnN(); b++) {
                    soe.HG[a][b] = 0;
                    soe.CG[a][b] = 0;
                    soe.PG[a] = 0;
                }
            }


            for (int l = 0; l < globalData.getnE(); l++) {

                temp_HBC = elem4.funHBC(l, nodeList, elementList, globalData.getSchematc(), globalData.getAlfa());
                temp_P = elem4.funP(l, nodeList, elementList, globalData.getSchematc(), globalData.getAlfa(), globalData.getTotoczenia());

                for (int i = 0; i < (int)Math.pow(globalData.getSchematc(),2); i++) {
                    temp_H = elem4.funH(i, l, nodeList, elementList, globalData.getK());
                    temp_C = elem4.funC(i, l, nodeList, elementList, globalData.getRo(), globalData.getC());

                    for (int j = 0; j < 4; j++) {
                        for (int k = 0; k < 4; k++) {
                            elementList.get(l).Hl[j][k] += (temp_H[j][k]);
                            elementList.get(l).Cl[j][k] += temp_C[j][k];
                        }
                    }
                }


                for (int j = 0; j < 4; j++) {
                    for (int k = 0; k < 4; k++) {
                        elementList.get(l).Hl[j][k] += temp_HBC[j][k];
                    }
                }

                /*System.out.println("Macierz HBC dla elementu: " + (l + 1));
                for (int j = 0; j < 4; j++) {
                    for (int k = 0; k < 4; k++) {
                        System.out.print(temp_HBC[j][k] + " ");
                    }
                    System.out.println();
                }
                System.out.println();*/

                /*System.out.println("Macierz Hl dla elementu: " + (l + 1));
                for (int j = 0; j < 4; j++) {
                    for (int k = 0; k < 4; k++) {
                        System.out.print(elementList.get(l).Hl[j][k] + " ");
                    }
                    System.out.println();
                }
                System.out.println();*/


                /*System.out.println("Macierz Cl dla elementu: " + (l + 1));
                for (int j = 0; j < 4; j++) {
                    for (int k = 0; k < 4; k++) {
                        System.out.print(elementList.get(l).Cl[j][k] + " ");
                    }
                    System.out.println();
                }
                System.out.println();*/

                //System.out.println("Pl dla elementu: " + (l + 1));
                for (int j = 0; j < 4; j++) {
                    elementList.get(l).Pl[j] = temp_P[j];
                    //System.out.print(elementList.get(l).Pl[j] + " ");
                }

            }


            //H globalne  i C globalne
            for (int l = 0; l < globalData.getnE(); l++) {
                for (int i = 0; i < 4; i++) {
                    for (int j = 0; j < 4; j++) {
                        soe.HG[elementList.get(l).id[i] - 1][elementList.get(l).id[j] - 1] += elementList.get(l).Hl[i][j];
                        soe.CG[elementList.get(l).id[i] - 1][elementList.get(l).id[j] - 1] += elementList.get(l).Cl[i][j];
                    }
                    soe.PG[elementList.get(l).id[i] - 1] += elementList.get(l).Pl[i];
                }
            }

            /*System.out.println("\nMacierz H globalna: ");
            for(int a = 0; a < globalData.getnN(); a++) {
                for (int b = 0; b < globalData.getnN(); b++) {
                    System.out.print(soe.HG[a][b] + ", ");
                }
                System.out.println();
            }

            System.out.println("\nMacierz C globalna: ");
            for(int a = 0; a < globalData.getnN(); a++){
                for(int b = 0; b < globalData.getnN(); b++){
                    System.out.print(soe.CG[a][b] + ", ");
                }
                System.out.println();
            }

            System.out.println("\nP globalne: ");
            for(int a = 0; a < globalData.getnN(); a++){
                    System.out.print(soe.PG[a] + ", ");
            }
            System.out.println();*/



            //System.out.println("\n[H]+[C]/dT ");

            for (int a = 0; a < globalData.getnN(); a++) {
                for (int b = 0; b < globalData.getnN(); b++) {
                    soe.HG[a][b] = (soe.HG[a][b] + (soe.CG[a][b] / globalData.getDeltat()));
                }
            }

            /*for (int a = 0; a < globalData.getnN(); a++) {
                for (int b = 0; b < globalData.getnN(); b++) {
                    System.out.print(soe.HG[a][b] + " ");
                }
                System.out.println();
            }*/


            //System.out.println("\n[P]+([C]/dT)*[T0]");

            for (int i = 0; i < globalData.getnN(); i++) {
                double s = 0;
                for (int j = 0; j < globalData.getnN(); j++) {
                    s += (soe.CG[i][j] / globalData.getDeltat()) * nodeList.get(j).getT0();
                }
                soe.PG[i] = (-1) * soe.PG[i] + s;
            }

            /*for (int a = 0; a < globalData.getnN(); a++) {
                System.out.print(soe.PG[a] + ", ");
            }
            System.out.println();*/


            soe.t1 = GaussElimination(soe.HG, soe.PG, globalData.getnN());

            /*System.out.println("\nt1: ");
            for (int a = 0; a < globalData.getnN(); a++) {
                System.out.println(soe.t1[a] + ", ");
            }*/

            double min = minimalny(soe.t1, globalData.getnN());
            double max = maksymalny(soe.t1, globalData.getnN());
            System.out.println("min: " + min + " max: " + max);
            System.out.println();

            for (int a = 0; a < globalData.getnN(); a++) {
                nodeList.get(a).setT0(soe.t1[a]);
            }
        }
    }

    public static double minimalny(double[] A, int n) {
        int i;
        double min;
        min = A[0];
        for (i = 0; i < n; i++)
            if (A[i] < min)
                min = A[i];
        return min;
    }

    public static double maksymalny(double[] A, int n) {
        int i;
        double max;
        max = A[0];
        for (i = 0; i < n; i++)
            if (A[i] > max)
                max = A[i];
        return max;
    }


    static public double[] GaussElimination(double[][] A, double[] b, int n) {
        double[] x = new double[n];

        double[][] tmpA = new double[n][n + 1];

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                tmpA[i][j] = A[i][j];
            }
            tmpA[i][n] = b[i];
        }
        double tmp = 0;
        for (int k = 0; k < n - 1; k++) {
            for (int i = k + 1; i < n; i++) {
                tmp = tmpA[i][k] / tmpA[k][k];
                for (int j = k; j < n + 1; j++) {
                    tmpA[i][j] -= tmp * tmpA[k][j];
                }
            }
        }
        for (int k = n - 1; k >= 0; k--) {
            tmp = 0;
            for (int j = k + 1; j < n; j++) {
                tmp += tmpA[k][j] * x[j];
            }
            x[k] = (tmpA[k][n] - tmp) / tmpA[k][k];
        }
        return x;
    }
}
