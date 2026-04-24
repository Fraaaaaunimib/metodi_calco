import org.ejml.data.DMatrixRMaj;
import org.ejml.data.DMatrixSparseCSC;
import org.ejml.dense.row.CommonOps_DDRM;
import org.ejml.dense.row.NormOps_DDRM;
import com.funzioni.Funzioni;
import java.io.File;
import java.util.Arrays;
import java.util.Scanner;
import java.io.IOException;
import java.text.DecimalFormat;


import org.ejml.sparse.csc.CommonOps_DSCC;

public class Jacobi {

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
        DMatrixRMaj diagValues = new DMatrixRMaj(A.numRows, 1);
        // Estrai i valori diagonali di A e copia i loro valori in diagValues
        CommonOps_DSCC.extractDiag(A, diagValues);
        // Crea matrice diagonale D con tutti 0 di tipo sparso
        DMatrixSparseCSC D = new DMatrixSparseCSC(A.numRows, A.numCols);
        // Riempie la matrice diagonale D con i valori diagonali di A
        CommonOps_DSCC.diag(D, diagValues.data, 0, A.numRows); 

        DMatrixRMaj dInv = new DMatrixRMaj(A.numRows, 1);
        for (int i = 0; i < A.numRows; i++) {
            dInv.set(i, 1.0 / diagValues.get(i, 0));
        }


        DMatrixSparseCSC B = new DMatrixSparseCSC(numRows, numCols);
        CommonOps_DSCC.add(1.0, D, -1.0, A, B, null, null); // Calcola B = D - A

        double[] xOld = x0.clone();
        double[] xNew = Arrays.copyOf(xOld, xOld.length+1);
        int nit = 0;


        long inizio = System.nanoTime();
        while(Funzioni.normaInfinito(xNew, xOld) > tol && nit < maxIter){
            xOld = xNew.clone();
            DMatrixRMaj bVec = new DMatrixRMaj(b);
            DMatrixRMaj temp = new DMatrixRMaj(numRows, 1);
            CommonOps_DDRM.mult(B, xOld, temp);
            CommonOps_DDRM.addEquals(temp, bVec);
            CommonOps_DDRM.elementMult(dInv, temp, temp);
            nit++;
        }
        long durata = System.nanoTime() - inizio;
        DecimalFormat secondi = new DecimalFormat().format(durata);

        DMatrixRMaj diff = new DMatrixRMaj(xNew.length, 1);
        double err = (NormOps_DDRM.normPInf(CommonOps_DDRM.subtract(xNew, xOld,diff)) / NormOps_DDRM.normPInf(xNew));

        return new double[]{Double.parseDouble(secondi), nit, err};
    }

}