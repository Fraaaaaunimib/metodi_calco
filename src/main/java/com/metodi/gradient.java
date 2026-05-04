package com.metodi;

import org.ejml.data.DMatrixRMaj;
import org.ejml.data.DMatrixSparseCSC;
import org.ejml.dense.row.CommonOps_DDRM;
import org.ejml.dense.row.NormOps_DDRM;
import org.ejml.sparse.csc.CommonOps_DSCC;

import com.funzioni.Risultato;

public class gradient{
    public static Risultato gradiente(DMatrixSparseCSC A, DMatrixRMaj b, DMatrixRMaj x0, int nMax, double tol){
        int n = A.numCols;
        int m = A.numRows;
        int L = x0.getNumRows();

        if (n != m) {
            System.out.println("La matrice A deve essere quadrata.");
            return null;
        } else if (L != m) {
            System.out.println("La dimensione di A deve essere uguale alla dimensione di x0.");
            return null;
        }

        DMatrixRMaj xK = x0.copy();
        DMatrixRMaj rK = new DMatrixRMaj(x0.numRows, 1);
        CommonOps_DSCC.mult(A, xK, rK);
        CommonOps_DDRM.subtract(b, rK, rK);
        double normaB = NormOps_DDRM.normP2(b);
        if(normaB == 0.0 || normaB == 1.0){
            return new Risultato(null, Integer.MAX_VALUE, Double.NaN, Double.NaN);
        }
        int nit = 0;
        double err = (NormOps_DDRM.normP2(rK)/normaB);
        Long inizio = System.nanoTime();

        while(err>tol && nit<nMax){
            DMatrixRMaj Ark = new DMatrixRMaj(A.numRows, 1);
            CommonOps_DSCC.mult(A, rK, Ark);
            DMatrixRMaj rKT = new DMatrixRMaj(1, x0.numRows);
            CommonOps_DDRM.transpose(rK, rKT);
            double step = (CommonOps_DDRM.dot(rK, rKT) / (CommonOps_DDRM.dot(Ark, rKT)));
            CommonOps_DDRM.add(xK, step, rK, xK);
            CommonOps_DDRM.add(rK, -step, Ark, rK);
            err = NormOps_DDRM.normP2(rK) / normaB;
            nit++;
        }

        Long durata = System.nanoTime() - inizio;
        double tempoSecondi = durata / 1_000_000_000.0;

        DMatrixRMaj xEsatta = new DMatrixRMaj(xK.getNumRows(), 1);
        for (int i = 0; i < xEsatta.getNumRows(); i++) {
            xEsatta.set(i, 1.0);
        }

        // Calcola differenza xNew - xEsatta
        DMatrixRMaj diffErr = new DMatrixRMaj(xK.getNumRows(), 1);
        CommonOps_DDRM.subtract(xK, xEsatta, diffErr);

        // Errore relativo
        err = NormOps_DDRM.normP2(diffErr) / NormOps_DDRM.normP2(xEsatta);

        return new Risultato(xK, nit, tempoSecondi, err);
    }
}