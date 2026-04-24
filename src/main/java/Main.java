import org.ejml.data.DMatrixRMaj;
import org.ejml.data.DMatrixSparseCSC;
import org.ejml.sparse.csc.CommonOps_DSCC;

import com.funzioni.*;
import com.metodi.Jacobi;

public class Main {
    public static void main(String[] args) {
        DMatrixSparseCSC matrice = Funzioni.caricaMatriceManualmente("spa1.mtx");
        double[] tolleranze = {1e-4, 1e-6, 1e-8, 1e-10};

        Jacobi jacobi = new Jacobi();


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

        for(double tol : tolleranze) {
            Risultato r = jacobi.jacobi(matrice, b, x0, 20000, tol);

            System.out.println("Tol = " + tol);
            System.out.println("Iterazioni = " + r.getNit());
            System.out.println("Errore = " + r.getErrore());
            System.out.println("Tempo = " + r.getTempo());
            System.out.println("----------------------");
        }
        

    }
}