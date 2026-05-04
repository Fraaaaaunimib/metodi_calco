package com.metodi;

import com.funzioni.*;

import org.ejml.data.DMatrixRMaj;
import org.ejml.data.DMatrixSparseCSC;
import org.ejml.dense.row.CommonOps_DDRM;
import org.ejml.dense.row.NormOps_DDRM;
import org.ejml.sparse.csc.CommonOps_DSCC;

public class Gauss {

    public static Risultato gauss(DMatrixSparseCSC A, DMatrixRMaj b, DMatrixRMaj x0, int maxIter, double tol) {

        int n = A.numCols;
        int m = A.numRows;

        if (m != n) {
            System.out.println("Matrice A non è quadrata");
            return null;
        }

        if (x0.numRows != n) {
            System.out.println("Dimensione di x0 errata");
            return null;
        }

        // =========================
        // PRE-CALCOLO
        // =========================
        DMatrixSparseCSC L = Funzioni.triang_inf(A); // L + D
        DMatrixSparseCSC U = new DMatrixSparseCSC(n, n, A.nz_length);

        // U = A - L
        CommonOps_DSCC.add(1.0, A, -1.0, L, U, null, null);

        // =========================
        // INIZIALIZZAZIONE
        // =========================
        DMatrixRMaj xNew = x0.copy();
        DMatrixRMaj xOld = new DMatrixRMaj(n, 1);

        DMatrixRMaj temp = new DMatrixRMaj(n, 1);
        DMatrixRMaj diff = new DMatrixRMaj(n, 1);
        DMatrixRMaj xDiff = new DMatrixRMaj(n, 1);

        int nit = 0;
        double normaDiff = Double.MAX_VALUE;

        long inizio = System.nanoTime();

        // =========================
        // CICLO ITERATIVO
        // =========================
        while (normaDiff > tol && nit < maxIter) {

            // xOld = xNew
            xOld = xNew.copy();

            // temp = U * xOld
            CommonOps_DSCC.mult(U, xOld, temp);

            // diff = b - temp
            CommonOps_DDRM.subtract(b, temp, diff);

            // risolvi L xNew = diff
            DMatrixRMaj xNext = Funzioni.solveTriangInf(L, diff);

            // calcolo differenza
            CommonOps_DDRM.subtract(xNext, xOld, xDiff);
            normaDiff = NormOps_DDRM.normPInf(xDiff);

            // aggiorna soluzione
            xNew = xNext.copy();

            nit++;
        }

        long durata = System.nanoTime() - inizio;
        double tempoSecondi = durata / 1_000_000_000.0;

        // =========================
        // ERRORE (se soluzione esatta = [1,...,1])
        // =========================
        DMatrixRMaj xEsatta = new DMatrixRMaj(n, 1);
        for (int i = 0; i < n; i++) {
            xEsatta.set(i, 0, 1.0);
        }

        DMatrixRMaj diffErr = new DMatrixRMaj(n, 1);
        CommonOps_DDRM.subtract(xNew, xEsatta, diffErr);

        double err = NormOps_DDRM.normPInf(diffErr) / NormOps_DDRM.normPInf(xEsatta);

        return new Risultato(xNew, nit, tempoSecondi, err);
    }
}