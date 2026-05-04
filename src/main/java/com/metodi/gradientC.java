package com.metodi;

import org.ejml.data.DMatrixRMaj;
import org.ejml.data.DMatrixSparseCSC;
import org.ejml.dense.row.CommonOps_DDRM;
import org.ejml.dense.row.NormOps_DDRM;
import org.ejml.sparse.csc.CommonOps_DSCC;

import com.funzioni.*;

public class gradientC{
    /*
    Implementazione del metodo del gradiente coniugato per matrici sparse
    Parametri: A (matrice sparsa di partenza), b (vettore di termini noti), x0 (vettore nullo), maxIter (numero massimo di iterazioni, deve essere >=20000), tol (tolleranza)
    Ritorno: oggetto di tipo Risultato (numero di iterazioni, tempo di esecuzione, errore relativo)
    Errore:
    - Se la matrice A non è quadrata
    - Se la lunghezza di x0 è diversa da A
    */
    public static Risultato gradienteC(DMatrixSparseCSC A, DMatrixRMaj b, DMatrixRMaj x0, int nMax, double tol){

        try {
            Funzioni.controlloDimensione(A.numCols, A.numRows, x0.numRows);
        } catch (Exception e){
            System.err.println(e);
        }

        // r = r^(k)
        DMatrixRMaj r = new DMatrixRMaj(b.numRows, 1);

        // r1 = r^(k+1)
        DMatrixRMaj r1 = new DMatrixRMaj(r);

        // r0 = b - A*x^(0)
        CommonOps_DSCC.mult(A, x0, r);
        CommonOps_DDRM.subtract(b, r, r);

        // p0 = p^(k)
        DMatrixRMaj p0 = r.copy();

        // pT = p^(k)_T
        DMatrixRMaj pT = new DMatrixRMaj(p0);
        
        int nit = 0;
        double err = (NormOps_DDRM.normP2(r)/NormOps_DDRM.normP2(b));

        Long inizio = System.nanoTime();

        DMatrixRMaj xK = x0.copy();

        while(err>tol && nit<nMax){
            // Apk = A*p^(k)
            DMatrixRMaj Apk = new DMatrixRMaj(A.numRows, 1);
            CommonOps_DSCC.mult(A, p0, Apk);

            CommonOps_DDRM.transpose(p0, pT);

            // step = ((r^(k)_T * r^(k)) / (p^(k)_T * A * p^(k)))
            double step = (CommonOps_DDRM.dot(r, r) / (CommonOps_DDRM.dot(Apk, pT)));

            // Aggiornamento della soluzione: x^(k+1) = x^(k) + step_k * p^(k)
            CommonOps_DDRM.add(xK, step, p0, xK);

            // Aggiornamento del residuo: r^(k+1) = r^(k) - step_k * A * p^(k)
            CommonOps_DDRM.add(r, -step, Apk, r1);

            err = NormOps_DDRM.normP2(r1) / NormOps_DDRM.normP2(b);
            nit++;

            // Coefficiente di correzione: beta_k = ((r^(k+1)_T * r^(k+1)) / (r^(k)_T * r^(k)))
            double beta = (CommonOps_DDRM.dot(r1, r1)/(CommonOps_DDRM.dot(r, r)));

            // Nuova direzione coniugata: p^(k+1) = r^(k+1) + beta_k * p^(k)
            CommonOps_DDRM.add(r1, beta, p0, p0);

            // Aggiornamento di r con il nuovo residuo per riuso al prossimo ciclo
            r = r1.copy();
        }

        Long durata = System.nanoTime() - inizio;
        double tempoSecondi = durata / 1_000_000_000.0;

        return new Risultato(nit, tempoSecondi, Funzioni.calcoloErroreRelativo(A.numRows, xK));

    }
}