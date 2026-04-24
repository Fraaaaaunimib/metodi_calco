package com.metodi;
import org.ejml.data.DMatrix;
import org.ejml.data.DMatrixRMaj;
import org.ejml.data.DMatrixSparseCSC;
import org.ejml.dense.row.CommonOps_DDRM;
import org.ejml.dense.row.NormOps_DDRM;
import com.funzioni.*;
import java.util.Arrays;

import org.ejml.sparse.csc.CommonOps_DSCC;

public class Jacobi {

    // Implementazione del metodo di Jacobi per risolvere il sistema lineare Ax = b
    public Risultato jacobi(final DMatrixSparseCSC A, final DMatrixRMaj b, final DMatrixRMaj x0, final int maxIter, final double tol) {
        int n = A.numCols;
        int m = A.numRows;
        int L = x0.getNumRows();
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

        DMatrixRMaj xOld = x0.copy();
        DMatrixRMaj xNew = xOld.copy();

        DMatrixRMaj diffOld = new DMatrixRMaj(xOld.getNumRows(), 1);
        Arrays.fill(diffOld.data, Double.MAX_VALUE);

       // CommonOps_DDRM.subtract(DMatrixRMaj.wrap(xNew.getNumRows(), 1, xNew.getData()), xOld, diffOld);

        int nit = 0;
        DMatrixRMaj bVec = new DMatrixRMaj(b);
        DMatrixRMaj temp = new DMatrixRMaj(A.numRows, 1);

        long inizio = System.nanoTime();
        while(NormOps_DDRM.normPInf(diffOld) > tol && nit < maxIter){
            xOld = xNew.copy();
            CommonOps_DSCC.mult(B, xOld, temp);
            CommonOps_DDRM.addEquals(temp, bVec);
            CommonOps_DDRM.elementMult(dInv, temp, temp);
            nit++;
            CommonOps_DDRM.subtract(temp, xOld, diffOld);
            xNew = temp.copy();
        }
        long durata = System.nanoTime() - inizio;
        double tempoSecondi = durata / 1_000_000_000.0;
        
        // Costruisci soluzione esatta: vettore di tutti 1
        DMatrixRMaj xEsatta = new DMatrixRMaj(xNew.getNumRows(), 1);
        for(int i = 0; i < xEsatta.getNumRows(); i++) {
            xEsatta.set(i, 1.0);
        }

        // Calcola differenza xNew - xEsatta
        DMatrixRMaj diffErr = new DMatrixRMaj(xNew.getNumRows(), 1);
        CommonOps_DDRM.subtract(xNew, xEsatta, diffErr);

        // Errore relativo
        double err = NormOps_DDRM.normPInf(diffErr) / NormOps_DDRM.normPInf(xEsatta);
        return new Risultato(xNew, nit, tempoSecondi, err);
    }

}