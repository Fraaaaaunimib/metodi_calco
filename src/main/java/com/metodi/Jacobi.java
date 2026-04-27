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
        for(int i = 0; i < diagValues.getNumRows(); i++) {
            dInv.set(i, 0, 1.0 / diagValues.get(i, 0));
        }

        DMatrixSparseCSC B = new DMatrixSparseCSC(A.numRows, A.numCols);
        CommonOps_DSCC.add(1.0, D, -1.0, A, B, null, null); // Calcola B = D - A

        DMatrixRMaj x = new DMatrixRMaj(x0);

        int nit = 0;
        DMatrixRMaj diff = new DMatrixRMaj(n, 1);
        DMatrixRMaj temp = new DMatrixRMaj(n, 1);

        double normaDiff = calcoloNormaDiff(A, x, b, diff);

        long inizio = System.nanoTime();
        
        while(normaDiff > tol && nit < maxIter){
            CommonOps_DSCC.mult(A, x, temp);
            CommonOps_DDRM.subtract(b, temp, temp);
            CommonOps_DDRM.elementDiv(temp, diagValues, temp);
            CommonOps_DDRM.addEquals(x, temp);

            nit++;
            normaDiff = calcoloNormaDiff(A, x, b, diff);
        }
        long durata = System.nanoTime() - inizio;
        double tempoSecondi = durata / 1_000_000_000.0;
        
        // Costruisci soluzione esatta: vettore di tutti 1
        DMatrixRMaj xEsatta = new DMatrixRMaj(x.getNumRows(), 1);
        for(int i = 0; i < xEsatta.getNumRows(); i++) {
            xEsatta.set(i, 1.0);
        }

        // Calcola differenza xNew - xEsatta
        DMatrixRMaj diffErr = new DMatrixRMaj(x.getNumRows(), 1);
        CommonOps_DDRM.subtract(x, xEsatta, diffErr);

        // Errore relativo
        double err = NormOps_DDRM.normPInf(diffErr) / NormOps_DDRM.normPInf(xEsatta);
        return new Risultato(x, nit, tempoSecondi, err);
    }

    public double calcoloNormaDiff(final DMatrixSparseCSC A, final DMatrixRMaj x0, final DMatrixRMaj b, final DMatrixRMaj workspace) {
        // ||A * x0 - b||_inf / ||b||_inf
        CommonOps_DSCC.mult(A, x0, workspace);
        CommonOps_DDRM.subtract(workspace, b, workspace);
        double normaDiff =NormOps_DDRM.normPInf(workspace);
        double normaInf = NormOps_DDRM.normPInf(b);
        normaDiff = normaDiff / normaInf;

        return normaDiff;
    }

}