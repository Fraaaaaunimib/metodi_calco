package com.metodi;

import org.ejml.data.DMatrixRMaj;
import org.ejml.data.DMatrixSparseCSC;
import org.ejml.dense.row.CommonOps_DDRM;
import com.funzioni.*;

import org.ejml.sparse.csc.CommonOps_DSCC;

public class Jacobi {

    /*
    Implementazione del metodo di Jacobi per matrici sparse
    Parametri: A (matrice sparsa di partenza), b (vettore di termini noti), x0 (vettore nullo), maxIter (numero massimo di iterazioni, deve essere >=20000), tol (tolleranza)
    Ritorno: oggetto di tipo Risultato (numero di iterazioni, tempo di esecuzione, errore relativo)
    Errore:
    - Se la matrice A non è quadrata
    - Se la lunghezza di x0 è diversa da A
    */
    public static Risultato jacobi(final DMatrixSparseCSC A, final DMatrixRMaj b, final DMatrixRMaj x0, final int maxIter,
            final double tol) {
        int n = A.numCols;
        int m = A.numRows;
        int L = x0.getNumRows();

        try {
            Funzioni.controlloDimensione(n, m, L);
            if (!Funzioni.isNotDiagonaleZero(A)) {
            throw new Exception("La matrice A deve avere elementi non zero sulla diagonale.");
            }
        } catch (Exception e){
            System.err.println(e);
        }

        // Crea array per i valori diagonali di A
        DMatrixRMaj diagValues = new DMatrixRMaj(A.numRows, 1);
        // Estrai i valori diagonali di A e copia i loro valori in diagValues
        CommonOps_DSCC.extractDiag(A, diagValues);
        // Crea matrice diagonale D con tutti 0 di tipo sparso
        DMatrixSparseCSC D = new DMatrixSparseCSC(A.numRows, A.numCols);
        // Riempie la matrice diagonale D con i valori diagonali di A
        CommonOps_DSCC.diag(D, diagValues.data, 0, A.numRows);

        // Calcolo dell'inversa della diagonale con 1/valore
        DMatrixRMaj dInv = new DMatrixRMaj(A.numRows, 1);
        for (int i = 0; i < diagValues.getNumRows(); i++) {
            dInv.set(i, 0, 1.0 / diagValues.get(i, 0));
        }

        // B = D - A
        DMatrixSparseCSC B = new DMatrixSparseCSC(A.numRows, A.numCols);
        CommonOps_DSCC.add(1.0, D, -1.0, A, B, null, null);

        DMatrixRMaj x = new DMatrixRMaj(x0);

        int nit = 0;
        DMatrixRMaj temp = new DMatrixRMaj(n, 1);

        double normaDiff = Funzioni.calcoloNormaDiff(A, x, b);

        long inizio = System.nanoTime();

        while (normaDiff > tol && nit < maxIter) {
            // x = ((b - A*x) / D) + x
            CommonOps_DSCC.mult(A, x, temp);
            CommonOps_DDRM.subtract(b, temp, temp);
            CommonOps_DDRM.elementDiv(temp, diagValues, temp);
            CommonOps_DDRM.addEquals(x, temp);
            
            nit++;

            normaDiff = Funzioni.calcoloNormaDiff(A, x, b);
        }
        long durata = System.nanoTime() - inizio;
        double tempoSecondi = durata / 1_000_000_000.0;

        return new Risultato(nit, tempoSecondi, Funzioni.calcoloErroreRelativo(n, x));
    }

}