package com.metodi;
import org.ejml.data.DMatrix;
import org.ejml.data.DMatrixRMaj;
import org.ejml.data.DMatrixSparseCSC;
import org.ejml.dense.row.CommonOps_DDRM;
import org.ejml.dense.row.NormOps_DDRM;
import com.funzioni.Funzioni;
import java.util.Arrays;

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


        DMatrixSparseCSC B = new DMatrixSparseCSC(A.numRows, A.numCols);
        CommonOps_DSCC.add(1.0, D, -1.0, A, B, null, null); // Calcola B = D - A

        DMatrixRMaj xOld = DMatrixRMaj.wrap(x0.length, 1, x0);
        DMatrixRMaj xNew = DMatrixRMaj.wrap(xOld.getNumRows(), 1, xOld.getData());
        for (int i = 0; i < xOld.getNumRows(); i++) {
            xNew.set(i, xOld.get(i)+1.0);
        }
        DMatrixRMaj diffOld = new DMatrixRMaj(xOld.getNumRows(), 1);

        CommonOps_DDRM.subtract(DMatrixRMaj.wrap(xNew.getNumRows(), 1, xNew.getData()), xOld, diffOld);

        int nit = 0;
        DMatrixRMaj bVec = new DMatrixRMaj(b);
        DMatrixRMaj temp = new DMatrixRMaj(A.numRows, 1);

        long inizio = System.nanoTime();
        while(NormOps_DDRM.normPInf(diffOld) > tol && nit < maxIter){
            for (int i = 0; i < xOld.getNumRows(); i++) {
                xNew.set(i, xOld.get(i));
            }
            CommonOps_DSCC.mult(B, xOld, temp);
            CommonOps_DDRM.addEquals(temp, bVec);
            CommonOps_DDRM.elementMult(dInv, temp, temp);
            nit++;
            CommonOps_DDRM.subtract(temp, xOld, diffOld); 
        }
        long durata = System.nanoTime() - inizio;
        double tempoSecondi = durata / 1_000_000_000.0;

        DMatrixRMaj diff = new DMatrixRMaj(xNew.getNumRows(), 1);
        double err = (NormOps_DDRM.normPInf(CommonOps_DDRM.subtract(DMatrixRMaj.wrap(xNew.getNumRows(), 1, xNew.getData()), DMatrixRMaj.wrap(xOld.getNumRows(), 1, xOld.getData()), diff)) / NormOps_DDRM.normPInf(DMatrixRMaj.wrap(xNew.getNumRows(), 1, xNew.getData())));

        return new double[]{tempoSecondi, nit, err};
    }

}