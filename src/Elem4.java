import java.util.List;

public class Elem4 {

    double[][] dNde;
    double[][] dNdn;
    double[][] tabN;

    double[] w1;
    double[] w2;

    double[] dxde;
    double[] dyde;
    double[] dydn;
    double[] dxdn;

    double[] e;
    double[] n;



    public Elem4(int schematc) {                      //inicjalizacja tablic pochodnych funkcji kształtu

        dNde = new double[(int)Math.pow(schematc,2)][4];
        dNdn = new double[(int)Math.pow(schematc,2)][4];
        tabN = new double[(int)Math.pow(schematc,2)][4];

        dxde = new double[(int)Math.pow(schematc,2)];
        dyde = new double[(int)Math.pow(schematc,2)];
        dydn = new double[(int)Math.pow(schematc,2)];
        dxdn = new double[(int)Math.pow(schematc,2)];

        if (schematc==2) {

            e = new double [] { -1.0 / Math.sqrt(3.0), 1.0 / Math.sqrt(3.0), 1.0 / Math.sqrt(3.0), -1.0 / Math.sqrt(3.0)};
            n = new double [] { -1.0 / Math.sqrt(3.0), -1.0 / Math.sqrt(3.0), 1.0 / Math.sqrt(3.0), 1.0 / Math.sqrt(3.0)};

            w1 = new double [] {1,1,1,1};
            w2 = new double [] {1,1,1,1};

            for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++) {
                    if (j == 0) {
                        dNde[i][j] = -0.25 * (1 - n[i]);
                        dNdn[i][j] = -0.25 * (1 - e[i]);
                    } else if (j == 1) {
                        dNde[i][j] = 0.25 * (1 - n[i]);
                        dNdn[i][j] = -0.25 * (1 + e[i]);
                    } else if (j == 2) {
                        dNde[i][j] = 0.25 * (1 + n[i]);
                        dNdn[i][j] = 0.25 * (1 + e[i]);
                    } else {
                        dNde[i][j] = -0.25 * (1 + n[i]);
                        dNdn[i][j] = 0.25 * (1 - e[i]);
                    }
                }
            }

            for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++) {
                    if (j == 0) tabN[i][j] = 0.25 *(1- e[i]) * (1 - n[i]);
                    if (j == 1) tabN[i][j] = 0.25 *(1+ e[i]) * (1 - n[i]);
                    if (j == 2) tabN[i][j] = 0.25 *(1+ e[i]) * (1 + n[i]);
                    if (j == 3) tabN[i][j] = 0.25 *(1- e[i]) * (1 + n[i]);
                }
            }

        }

        else if(schematc==3){

            e = new double [] { -1*Math.sqrt(3.0/5.0), 0, Math.sqrt(3.0/5.0), -1*Math.sqrt(3.0/5.0), 0, Math.sqrt(3.0/5.0), -1*Math.sqrt(3.0/5.0), 0, Math.sqrt(3.0/5.0)};
            n = new double [] { -1*Math.sqrt(3.0/5.0), -1*Math.sqrt(3.0/5.0), -1*Math.sqrt(3.0/5.0), 0,0,0, Math.sqrt(3.0/5.0), Math.sqrt(3.0/5.0), Math.sqrt(3.0/5.0)};


            w1 = new double [] { 5.0/9.0, 8.0/9.0, 5.0/9.0, 5.0/9.0, 8.0/9.0, 5.0/9.0, 5.0/9.0, 8.0/9.0, 5.0/9.0 };
            w2 = new double [] { 5.0/9.0, 5.0/9.0, 5.0/9.0, 8.0/9.0, 8.0/9.0, 8.0/9.0, 5.0/9.0, 5.0/9.0, 5.0/9.0 };


            for (int i = 0; i < 9; i++) {
                for (int j = 0; j < 4; j++) {
                    if (j == 0) {
                        dNde[i][j] = -0.25 * (1 - n[i]);
                        dNdn[i][j] = -0.25 * (1 - e[i]);
                    } else if (j == 1) {
                        dNde[i][j] = 0.25 * (1 - n[i]);
                        dNdn[i][j] = -0.25 * (1 + e[i]);
                    } else if (j == 2) {
                        dNde[i][j] = 0.25 * (1 + n[i]);
                        dNdn[i][j] = 0.25 * (1 + e[i]);
                    } else {
                        dNde[i][j] = -0.25 * (1 + n[i]);
                        dNdn[i][j] = 0.25 * (1 - e[i]);
                    }
                }
            }

            for (int i = 0; i < 9; i++) {
                for (int j = 0; j < 4; j++) {
                    if (j == 0) tabN[i][j] = 0.25 *(1- e[i]) * (1 - n[i]);
                    if (j == 1) tabN[i][j] = 0.25 *(1+ e[i]) * (1 - n[i]);
                    if (j == 2) tabN[i][j] = 0.25 *(1+ e[i]) * (1 + n[i]);
                    if (j == 3) tabN[i][j] = 0.25 *(1- e[i]) * (1 + n[i]);
                }
            }

        }
        else if(schematc==4){

            e = new double [] {
                -(Math.sqrt((3.0 / 7.0) + ((2.0 / 7.0) * Math.sqrt(6.0 / 5.0)))), -(Math.sqrt((3.0 / 7.0) - ((2.0 / 7.0) * Math.sqrt(6.0 / 5.0)))), Math.sqrt((3.0 / 7.0) - ((2.0 / 7.0) * Math.sqrt(6.0 / 5.0))), Math.sqrt((3.0 / 7.0) + ((2.0 / 7.0) * Math.sqrt(6.0 / 5.0))),
                -(Math.sqrt((3.0 / 7.0) + ((2.0 / 7.0) * Math.sqrt(6.0 / 5.0)))), -(Math.sqrt((3.0 / 7.0) - ((2.0 / 7.0) * Math.sqrt(6.0 / 5.0)))),Math.sqrt((3.0 / 7.0) - ((2.0 / 7.0) * Math.sqrt(6.0 / 5.0))), Math.sqrt((3.0 / 7.0) + ((2.0 / 7.0) * Math.sqrt(6.0 / 5.0))),
                -(Math.sqrt((3.0 / 7.0) + ((2.0 / 7.0) * Math.sqrt(6.0 / 5.0)))), -(Math.sqrt((3.0 / 7.0) - ((2.0 / 7.0) * Math.sqrt(6.0 / 5.0)))), Math.sqrt((3.0 / 7.0) - ((2.0 / 7.0) * Math.sqrt(6.0 / 5.0))), Math.sqrt((3.0 / 7.0) + ((2.0 / 7.0) * Math.sqrt(6.0 / 5.0))),
                -(Math.sqrt((3.0 / 7.0) + ((2.0 / 7.0) * Math.sqrt(6.0 / 5.0)))), -(Math.sqrt((3.0 / 7.0) - ((2.0 / 7.0) * Math.sqrt(6.0 / 5.0)))), Math.sqrt((3.0 / 7.0) - ((2.0 / 7.0) * Math.sqrt(6.0 / 5.0))), Math.sqrt((3.0 / 7.0) + ((2.0 / 7.0) * Math.sqrt(6.0 / 5.0)))};
            n = new double [] {
                    -(Math.sqrt((3.0 / 7.0) + ((2.0 / 7.0) * Math.sqrt(6.0 / 5.0)))), -(Math.sqrt((3.0 / 7.0) + ((2.0 / 7.0) * Math.sqrt(6.0 / 5.0)))), -(Math.sqrt((3.0 / 7.0) + ((2.0 / 7.0) * Math.sqrt(6.0 / 5.0)))), -(Math.sqrt((3.0 / 7.0) + ((2.0 / 7.0) * Math.sqrt(6.0 / 5.0)))),
                    -(Math.sqrt((3.0 / 7.0) - ((2.0 / 7.0) * Math.sqrt(6.0 / 5.0)))), -(Math.sqrt((3.0 / 7.0) - ((2.0 / 7.0) * Math.sqrt(6.0 / 5.0)))),-(Math.sqrt((3.0 / 7.0) - ((2.0 / 7.0) * Math.sqrt(6.0 / 5.0)))), -(Math.sqrt((3.0 / 7.0) - ((2.0 / 7.0) * Math.sqrt(6.0 / 5.0)))),
                    Math.sqrt((3.0 / 7.0) - ((2.0 / 7.0) * Math.sqrt(6.0 / 5.0))), Math.sqrt((3.0 / 7.0) - ((2.0 / 7.0) * Math.sqrt(6.0 / 5.0))), Math.sqrt((3.0 / 7.0) - ((2.0 / 7.0) * Math.sqrt(6.0 / 5.0))), Math.sqrt((3.0 / 7.0) - ((2.0 / 7.0) * Math.sqrt(6.0 / 5.0))),
                    Math.sqrt((3.0 / 7.0) + ((2.0 / 7.0) * Math.sqrt(6.0 / 5.0))), Math.sqrt((3.0 / 7.0) + ((2.0 / 7.0) * Math.sqrt(6.0 / 5.0))), Math.sqrt((3.0 / 7.0) + ((2.0 / 7.0) * Math.sqrt(6.0 / 5.0))), Math.sqrt((3.0 / 7.0) + ((2.0 / 7.0) * Math.sqrt(6.0 / 5.0)))};

            w1 = new double[]{(18.0 - Math.sqrt(30.0)) / 36.0, (18.0 + Math.sqrt(30.0)) / 36.0, (18.0 + Math.sqrt(30.0)) / 36.0, (18.0 - Math.sqrt(30.0)) / 36.0, (18.0 - Math.sqrt(30.0)) / 36.0, (18.0 + Math.sqrt(30.0)) / 36.0,
                    (18.0 + Math.sqrt(30.0)) / 36.0, (18.0 - Math.sqrt(30.0)) / 36.0, (18.0 - Math.sqrt(30.0)) / 36.0, (18.0 + Math.sqrt(30.0)) / 36.0, (18.0 + Math.sqrt(30.0)) / 36.0, (18.0 - Math.sqrt(30.0)) / 36.0,
                    (18.0 - Math.sqrt(30.0)) / 36.0, (18.0 + Math.sqrt(30.0)) / 36.0, (18.0 + Math.sqrt(30.0)) / 36.0, (18.0 - Math.sqrt(30.0)) / 36.0};
            w2 = new double[]{(18.0 - Math.sqrt(30.0)) / 36.0, (18.0 - Math.sqrt(30.0)) / 36.0, (18.0 - Math.sqrt(30.0)) / 36.0, (18.0 - Math.sqrt(30.0)) / 36.0, (18.0 + Math.sqrt(30.0)) / 36.0,
                    (18.0 + Math.sqrt(30.0)) / 36.0, (18.0 + Math.sqrt(30.0)) / 36.0, (18.0 + Math.sqrt(30.0)) / 36.0, (18.0 + Math.sqrt(30.0)) / 36.0, (18.0 + Math.sqrt(30.0)) / 36.0,
                    (18.0 + Math.sqrt(30.0)) / 36.0, (18.0 + Math.sqrt(30.0)) / 36.0, (18.0 - Math.sqrt(30.0)) / 36.0, (18.0 - Math.sqrt(30.0)) / 36.0, (18.0 - Math.sqrt(30.0)) / 36.0, (18.0 - Math.sqrt(30.0)) / 36.0};

            for (int i = 0; i < 16; i++) {
                for (int j = 0; j < 4; j++) {
                    if (j == 0) {
                        dNde[i][j] = -0.25 * (1 - n[i]);
                        dNdn[i][j] = -0.25 * (1 - e[i]);
                    } else if (j == 1) {
                        dNde[i][j] = 0.25 * (1 - n[i]);
                        dNdn[i][j] = -0.25 * (1 + e[i]);
                    } else if (j == 2) {
                        dNde[i][j] = 0.25 * (1 + n[i]);
                        dNdn[i][j] = 0.25 * (1 + e[i]);
                    } else {
                        dNde[i][j] = -0.25 * (1 + n[i]);
                        dNdn[i][j] = 0.25 * (1 - e[i]);
                    }
                }
            }

            for (int i = 0; i < 16; i++) {
                for (int j = 0; j < 4; j++) {
                    if (j == 0) tabN[i][j] = 0.25 *(1- e[i]) * (1 - n[i]);
                    if (j == 1) tabN[i][j] = 0.25 *(1+ e[i]) * (1 - n[i]);
                    if (j == 2) tabN[i][j] = 0.25 *(1+ e[i]) * (1 + n[i]);
                    if (j == 3) tabN[i][j] = 0.25 *(1- e[i]) * (1 + n[i]);
                }
            }

        }
    }

    public double[][] jacobiMatrix(int i, int e, List<Node> nodeList, List<Element> elementList){

        double [][] jMatrix = new double[2][2];

        dxde[i] = 0;
        dydn[i] = 0;
        dxdn[i] = 0;
        dyde[i] = 0;

        for(int j = 0; j < 4; j++) {
            dxde[i] = dxde[i] + dNde[i][j] * nodeList.get(elementList.get(e).id[j] - 1).getX();
            dydn[i] = dydn[i] + dNdn[i][j] * nodeList.get(elementList.get(e).id[j] - 1).getY();
            dxdn[i] = dxdn[i] + dNdn[i][j] * nodeList.get(elementList.get(e).id[j] - 1).getX();
            dyde[i] = dyde[i] + dNde[i][j] * nodeList.get(elementList.get(e).id[j] - 1).getY();
        }

        for(int j = 0; j<2;j++) {
            for(int k = 0; k<2;k++) {
                if(j==0 && k==0) jMatrix[j][k] = dxde[i];
                if(j==0 && k==1) jMatrix[j][k] = dyde[i];
                if(j==1 && k==0) jMatrix[j][k] = dxdn[i];
                if(j==1 && k==1) jMatrix[j][k] = dydn[i];
            }
        }

        return jMatrix;
    }

    public double[][] funH(int i, int e, List<Node> nodeList, List<Element> elementList, double kk) {

        double[] dNdx = new double[4];
        double[] dNdy = new double[4];

        double [][] jacobiMatrix = new double[2][2];
        double detJacobiMatrix=0;
        double [][] reverseJacobiMatrix = new double[2][2];

        jacobiMatrix = jacobiMatrix(i, e, nodeList, elementList);
        detJacobiMatrix = jacobiMatrix[0][0] * jacobiMatrix[1][1] - jacobiMatrix[0][1] * jacobiMatrix[1][0];

        for(int j = 0; j<2;j++){
            for(int k = 0; k<2;k++){
                if(j==0 && k==0) reverseJacobiMatrix[j][k] = (jacobiMatrix[1][1]) / detJacobiMatrix;
                if(j==0 && k==1) reverseJacobiMatrix[j][k] = (-jacobiMatrix[j][k]) / detJacobiMatrix;
                if(j==1 && k==0) reverseJacobiMatrix[j][k] = (-jacobiMatrix[j][k]) / detJacobiMatrix;
                if(j==1 && k==1) reverseJacobiMatrix[j][k] = (jacobiMatrix[0][0]) / detJacobiMatrix;
            }
        }

        for (int j = 0; j < 4; j++) {
            dNdx[j] = reverseJacobiMatrix[0][0] * dNde[i][j] + reverseJacobiMatrix[0][1] * dNdn[i][j];
            dNdy[j] = reverseJacobiMatrix[1][0] * dNde[i][j] + reverseJacobiMatrix[1][1] * dNdn[i][j];
        }

        double[][] dNdxdNdxT = new double[4][4];
        double[][] dNdydNdyT = new double[4][4];

        for (int m = 0; m < 4; m++) {
            for (int n = 0; n < 4; n++) {
                dNdxdNdxT[m][n] = dNdx[m] * dNdx[n];
                dNdydNdyT[m][n] = dNdy[m] * dNdy[n];
            }
        }

        double[][] H = new double[4][4];

        for (int m = 0; m < 4; m++) {
            for (int n = 0; n < 4; n++) {
                H[m][n] = (dNdxdNdxT[m][n] + dNdydNdyT[m][n]) * kk * detJacobiMatrix * w1[i] * w2[i];
            }
        }
        return H;
    }

    public double[][] funHBC( int el, List<Node> nodeList, List<Element> elementList, int lpktcal, double alfa){

        double [][] HBC = new double[4][4];
        double [] Np = new double [4];
        for(int i=0; i<4; i++){
            Node nodeA = nodeList.get(elementList.get(el).id[i] - 1);
            Node nodeB = nodeList.get(elementList.get(el).id[(i+1)%4] - 1);

            if(!(nodeA.getBC()==1 && nodeB.getBC()==1)) continue;

            for(int j=0; j<lpktcal; j++) {
                double ee = 0;
                double nn = 0;

                if (i==0){ ee = e[j]; nn = -1.0; }
                if (i==1){ ee = 1.0; nn = e[j]; }
                if (i==2){ ee = (-1.0)*e[j];nn = 1.0; }
                if (i==3){ ee = -1.0;nn =(-1.0) * e[j]; }

                for (int k = 0; k < 4; k++) {
                    if (k == 0) Np[k] = 0.25 *(1- ee) * (1 - nn);
                    if (k == 1) Np[k] = 0.25 *(1+ ee) * (1 - nn);
                    if (k == 2) Np[k] = 0.25 *(1+ ee) * (1 + nn);
                    if (k == 3) Np[k] = 0.25 *(1- ee) * (1 + nn);
                }
                double detJ = (Math.sqrt((nodeA.getX()-nodeB.getX())*(nodeA.getX()-nodeB.getX()) + (nodeA.getY()-nodeB.getY())*(nodeA.getY()-nodeB.getY())))/2;
                for(int k = 0; k<4;k++) {
                    for(int l = 0; l<4;l++){
                        HBC[k][l] += Np[k] * Np[l] * w1[j] * alfa * detJ;
                    }
                }
            }
        }

        return HBC;
    }

    public double[] funP( int el, List<Node> nodeList, List<Element> elementList, int lpktcal, double alfa, double totoczenia){

        double [] P = new double[4];
        double [] Np = new double [4];

        for(int i=0; i<4; i++){
            Node nodeA = nodeList.get(elementList.get(el).id[i] - 1);
            Node nodeB = nodeList.get(elementList.get(el).id[(i+1)%4] - 1);

            if(!(nodeA.getBC()==1 && nodeB.getBC()==1)) continue;
            for(int j=0; j<lpktcal; j++) {
                double ee = 0;        //e
                double nn = 0;

                if (i==0){ ee = e[j]; nn = -1.0; }
                if (i==1){ ee = 1.0; nn = e[j]; }
                if (i==2){ ee = (-1.0)*e[j];nn = 1.0; }
                if (i==3){ ee = -1.0;nn =(-1.0) * e[j]; }
                for (int k = 0; k < 4; k++) {
                    if (k == 0) Np[k] = 0.25 *(1- ee) * (1 - nn);
                    if (k == 1) Np[k] = 0.25 *(1+ ee) * (1 - nn);
                    if (k == 2) Np[k] = 0.25 *(1+ ee) * (1 + nn);
                    if (k == 3) Np[k] = 0.25 *(1- ee) * (1 + nn);
                }
                double detJ = (Math.sqrt((nodeA.getX()-nodeB.getX())*(nodeA.getX()-nodeB.getX()) + (nodeA.getY()-nodeB.getY())*(nodeA.getY()-nodeB.getY())))/2;
                for(int k = 0; k<4;k++) {
                        P[k] += ( (-1)*Np[k] * w1[j] * totoczenia * alfa * detJ);     //
                }
            }
        }
        return P;
    }



    public double[][] funC(int i, int e, List<Node> nodeList, List<Element> elementList, double ro, double c){

        double [][] jacobiMatrix = new double[2][2];
        double detJacobiMatrix=0;
        jacobiMatrix = jacobiMatrix(i, e, nodeList, elementList);
        detJacobiMatrix = jacobiMatrix[0][0] * jacobiMatrix[1][1] - jacobiMatrix[0][1] * jacobiMatrix[1][0];

        double [][] N_NT = new double[4][4];

        for(int k = 0; k<4;k++) {
            for(int j = 0; j<4;j++){
                N_NT[k][j] = tabN[i][k]*tabN[i][j];
            }
        }

        double [][] C = new double[4][4];

        for(int k = 0; k<4;k++){
            for(int j = 0; j<4; j++){
                C[k][j] = ((N_NT[k][j]) * c * ro * detJacobiMatrix * w1[i] * w2[i]);
            }
        }

        return C;
    }
}