package com.metodi;
import com.funzioni.*;

import java.util.Arrays;

import org.ejml.data.DMatrix;
import org.ejml.data.DMatrixRMaj;
import org.ejml.data.DMatrixSparseCSC;
import org.ejml.dense.row.CommonOps_DDRM;
import org.ejml.dense.row.NormOps_DDRM;
import org.ejml.sparse.csc.CommonOps_DSCC;

public class Gauss{
    public static Risultato gauss(DMatrixSparseCSC A, DMatrixRMaj b, DMatrixRMaj x0, int maxIter, double tol) {
    int n = A.numCols;
    int m = A.numRows;
    int L = x0.numRows;

    if(m != n){
        System.out.println("Matrice A non è quadrta");
        return null;
    }
    else if(L!=m){
        System.out.println("Dimensione di A diversa da x0");
        return null;
    }

    // 1. PRE-CALCOLO (Fuori dal while per la massima velocità)
    DMatrixSparseCSC A_triang = new DMatrixSparseCSC(n, n);
    A_triang = Funzioni.triang_inf(A);
    
    DMatrixSparseCSC B = new DMatrixSparseCSC(n, n);
    CommonOps_DSCC.add(1.0, A, -1.0, A_triang, B, null, null);

    // 2. INIZIALIZZAZIONE DATI
    DMatrixRMaj xNew = x0.copy();
    DMatrixRMaj xOld = x0.copy();
    DMatrixRMaj diff = new DMatrixRMaj(n, 1);
    diff.fill(0.0);
    
    // xnew = xold + 1 (come da schema prof)
    for(int i=0; i<n; i++){
        xNew.set(i, xOld.get(i) + 1.0);
    }

    int nit = 0;
    long inizio = System.nanoTime();
    // Prima del while, inizializza temp una volta sola
    DMatrixRMaj temp = new DMatrixRMaj(b.numRows, 1);

    double normaDiff = Funzioni.calcoloNormaDiff(A, xOld, b, diff);

    // 3. LOOP ITERATIVO
    while (normaDiff > tol && nit < maxIter) {
        xOld = xNew.copy();

        CommonOps_DSCC.mult(B, xOld, temp);
        CommonOps_DDRM.subtract(b, temp, diff);
        xNew = Funzioni.solveTriangInf(A_triang, diff);
        normaDiff = Funzioni.calcoloNormaDiff(A, xNew, b, diff);
        nit++;
    }

     long durata = System.nanoTime() - inizio;
        double tempoSecondi = durata / 1_000_000_000.0;

    DMatrixRMaj xEsatta = new DMatrixRMaj(xOld.getNumRows(), 1);
        for (int i = 0; i < xEsatta.getNumRows(); i++) {
            xEsatta.set(i, 1.0);
        }

        // Calcola differenza xNew - xEsatta
        DMatrixRMaj diffErr = new DMatrixRMaj(xOld.getNumRows(), 1);
        CommonOps_DDRM.subtract(xNew, xEsatta, diffErr);

        // Errore relativo
        double err = NormOps_DDRM.normP2(diffErr) / NormOps_DDRM.normP2(xEsatta);
        return new Risultato(xNew, nit, tempoSecondi, err);
}
}