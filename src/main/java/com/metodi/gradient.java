package com.metodi;

import org.ejml.data.DMatrixRMaj;
import org.ejml.data.DMatrixSparseCSC;
import org.ejml.dense.row.CommonOps_DDRM;
import org.ejml.dense.row.NormOps_DDRM;
import org.ejml.sparse.csc.CommonOps_DSCC;

import com.funzioni.*;

public class gradient{

    /*
    Implementazione del metodo del gradiente per matrici sparse
    Parametri: A (matrice sparsa di partenza), b (vettore di termini noti), x0 (vettore nullo), maxIter (numero massimo di iterazioni, deve essere >=20000), tol (tolleranza)
    Ritorno: oggetto di tipo Risultato (numero di iterazioni, tempo di esecuzione, errore relativo)
    Errore:
    - Se la matrice A non è quadrata
    - Se la lunghezza di x0 è diversa da A
    */
    public static Risultato gradiente(DMatrixSparseCSC A, DMatrixRMaj b, DMatrixRMaj x0, int nMax, double tol){
        int n = A.numCols;
        int m = A.numRows;
        int L = x0.getNumRows();

        try {
            Funzioni.controlloDimensione(n, m, L);
        } catch (Exception e){
            System.err.println(e);
        }

        DMatrixRMaj xK = x0.copy();
        DMatrixRMaj rK = new DMatrixRMaj(x0.numRows, 1);
        // Calcolo direzione: rK = b - A*xK 
        CommonOps_DSCC.mult(A, xK, rK);
        CommonOps_DDRM.subtract(b, rK, rK);

        int nit = 0;
        
        double err = (NormOps_DDRM.normP2(rK)/NormOps_DDRM.normP2(b));
        
        Long inizio = System.nanoTime();

        // Calcolo del passo step e aggiornamento della soluzione xK
        while(err>tol && nit<nMax){
            // Ark = A*r^(k)
            DMatrixRMaj Ark = new DMatrixRMaj(A.numRows, 1);
            CommonOps_DSCC.mult(A, rK, Ark);

            // rKT = r^(k)_T, trasposta di r^(k)
            DMatrixRMaj rKT = new DMatrixRMaj(1, x0.numRows);
            CommonOps_DDRM.transpose(rK, rKT);

            // step = ((r^(k)_T * r(k)) / (r^(k)_T * A*r^(k)))
            double step = (CommonOps_DDRM.dot(rK, rKT) / (CommonOps_DDRM.dot(Ark, rKT)));

            // Aggiorna la soluzione attuale x^(k+1) = x^(k) + step*r^(k)
            CommonOps_DDRM.add(xK, step, rK, xK);
            // Aggiorna il residuo attuale r^(k+1) = r^(k) - step*A*r^(k)
            CommonOps_DDRM.add(rK, -step, Ark, rK);
            err = NormOps_DDRM.normP2(rK) / NormOps_DDRM.normP2(b);
            nit++;
        }

        Long durata = System.nanoTime() - inizio;
        double tempoSecondi = durata / 1_000_000_000.0;

        return new Risultato(nit, tempoSecondi, Funzioni.calcoloErroreRelativo(n, xK));
    }
}