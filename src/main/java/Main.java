import org.ejml.data.DMatrixRMaj;
import org.ejml.data.DMatrixSparseCSC;
import org.ejml.sparse.csc.CommonOps_DSCC;

import com.funzioni.*;
import com.metodi.Jacobi;

public class Main {
    public static void main(String[] args) {
        DMatrixSparseCSC spa1 = Funzioni.caricaMatriceManualmente("matrici/spa1.mtx");
        DMatrixSparseCSC spa2 = Funzioni.caricaMatriceManualmente("matrici/spa2.mtx");
        DMatrixSparseCSC vem1 = Funzioni.caricaMatriceManualmente("matrici/vem1.mtx");
        DMatrixSparseCSC vem2 = Funzioni.caricaMatriceManualmente("matrici/vem2.mtx");

        String[] nomiMatrici = {"spa1", "spa2", "vem1", "vem2"};

        double[] tolleranze = {1e-4, 1e-6, 1e-8, 1e-10};

        Jacobi jacobi = new Jacobi();

        int numMatrice = 0;

        for (DMatrixSparseCSC matrice : new DMatrixSparseCSC[]{spa1, spa2, vem1, vem2}) {
            DMatrixRMaj x0 = new DMatrixRMaj(matrice.numRows, 1);
            for(int i = 0; i < x0.getNumRows(); i++) {
                x0.set(i, 0.0);
            }

            DMatrixRMaj x = new DMatrixRMaj(matrice.numRows, 1);
            for(int i = 0; i < x.getNumRows(); i++) {
                x.set(i, 1.0);
            }

            DMatrixRMaj b = new DMatrixRMaj(matrice.numRows, 1);
            CommonOps_DSCC.mult(matrice, x, b);

            System.out.println("Matrice: " + nomiMatrici[numMatrice]);
            for(double tol : tolleranze) {
            Risultato r = jacobi.jacobi(matrice, b, x0, 10000, tol);

            System.out.print("| Tol = " + tol);
            System.out.print(" | Iterazioni = " + r.getNit());
            System.out.print(" | Errore = " + r.getErrore());
            System.out.println(" | Tempo = " + r.getTempo() + " | ");
            System.out.println("----------------------");
            }
            System.out.println("");
            numMatrice++;
    }

    }
}