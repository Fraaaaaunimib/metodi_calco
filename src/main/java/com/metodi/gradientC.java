package com.metodi;

import org.ejml.data.DMatrixRMaj;
import org.ejml.data.DMatrixSparseCSC;
import org.ejml.dense.row.CommonOps_DDRM;
import org.ejml.dense.row.NormOps_DDRM;
import org.ejml.sparse.csc.CommonOps_DSCC;

import com.funzioni.Risultato;

public class gradientC{
    public static Risultato gradienteC(DMatrixSparseCSC A, DMatrixRMaj b, DMatrixRMaj x0, int nMax, double tol){
        DMatrixRMaj r = new DMatrixRMaj(b.numRows, 1);
        DMatrixRMaj r1 = new DMatrixRMaj(r);
        CommonOps_DSCC.mult(A, x0, r);
        CommonOps_DDRM.subtract(b, r, r);
        DMatrixRMaj p0 = r.copy();
        DMatrixRMaj pT = new DMatrixRMaj(p0);
        double normaB = NormOps_DDRM.normP2(b);
        if(normaB == 0.0 || normaB == 1.0){
            return new Risultato(null, Integer.MAX_VALUE, Double.NaN, Double.NaN);
        }
        int nit = 0;
        double err = (NormOps_DDRM.normP2(r)/normaB);
        Long inizio = System.nanoTime();
        DMatrixRMaj xK = x0.copy();

        while(err>tol && nit<nMax){
            DMatrixRMaj Apk = new DMatrixRMaj(A.numRows, 1);
            CommonOps_DSCC.mult(A, p0, Apk);
            CommonOps_DDRM.transpose(p0, pT);
            double step = (CommonOps_DDRM.dot(r, r) / (CommonOps_DDRM.dot(Apk, pT)));
            CommonOps_DDRM.add(xK, step, p0, xK);
            CommonOps_DDRM.add(r, -step, Apk, r1);
            err = NormOps_DDRM.normP2(r1) / normaB;
            nit++;
            double beta = (CommonOps_DDRM.dot(r1, r1)/(CommonOps_DDRM.dot(r, r)));
            CommonOps_DDRM.add(r1, beta, p0, p0);
            r = r1.copy();
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