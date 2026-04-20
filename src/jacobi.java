import org.ejml.data.DMatrixRMaj;
import org.ejml.data.DMatrixSparseCSC;
import java.io.File;
import java.util.Scanner;
import java.io.IOException;
import org.ejml.sparse.csc.CommonOps_DSCC;

public class jacobi {
    DMatrixSparseCSC matrice = Main.caricaMatriceManualmente("spa1.mtx");
    // Implementazione del metodo di Jacobi per risolvere il sistema lineare Ax = b
    public double[] jacobi(DMatrixSparseCSC A, double[] b, double[] x0, int maxIter, double tol) {
        int n = A.numCols;
        int m = A.numRows;
        int L = x0.length;
        if(n!=m){
            System.out.println("La matrice A deve essere quadrata.");
            return null;
        }
        else if(L!=m){
            System.out.println("La dimensione di A deve essere uguale alla dimensione di x0.");
            return null;
        }
        
        if(!Funzioni.diagonaleZeri(A)){
            System.out.println("La matrice A deve avere elementi non zero sulla diagonale.");
            return null;
        }

        // Crea array per i valori diagonali di A
        DMatrixRMaj diagValues = new DMatrixRMaj(numRows, 1);
        // Estrai i valori diagonali di A e copia i loro valori in diagValues
        CommonOps_DSCC.extractDiag(A, diagValues);
        // Crea matrice diagonale D con tutti 0 di tipo sparso
        DMatrixSparseCSC D = new DMatrixSparseCSC(A.numRows, A.numCols);
        // Riempie la matrice diagonale D con i valori diagonali di A
        CommonOps_DSCC.diag(D, diagValues.data, 0, A.numRows); 

        DMatrixSparseCSC B = new DMatrixSparseCSC(numRows, numCols);
        CommonOps_DSCC.add(1.0, D, -1.0, A, B, null, null); // Calcola B = D - A

    }

}