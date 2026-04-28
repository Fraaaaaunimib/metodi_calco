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

    // 1. PRE-CALCOLO (Fuori dal while per la massima velocità)
    DMatrixSparseCSC L = new DMatrixSparseCSC(n, n);
    L = Funzioni.triang_inf(A);
    
    DMatrixSparseCSC B = new DMatrixSparseCSC(n, n);
    CommonOps_DSCC.add(1.0, A, -1.0, L, B, null, null);

    // Salviamo le diagonali in un array per non cercarle più
    double[] diags = new double[n];
    for (int j = 0; j < n; j++) {
        diags[j] = L.get(j, j);
        if (Math.abs(diags[j]) < 1e-16) diags[j] = 1e-16; // Evita divisione per zero secca
    }

    // 2. INIZIALIZZAZIONE DATI
    DMatrixRMaj xNew = x0.copy();
    DMatrixRMaj xOld = new DMatrixRMaj(n, 1);
    DMatrixRMaj diff = new DMatrixRMaj(n, 1);
    
    // xnew = xold + 1 (come da schema prof)
    for(int i=0; i<n; i++) xOld.set(i, xNew.get(i) + 1.0);

    int nit = 0;
    long inizio = System.nanoTime();
       // Prima del while, inizializza temp una volta sola
DMatrixRMaj temp = new DMatrixRMaj(b.numRows, 1);
    // 3. LOOP ITERATIVO
    while (nit < maxIter) {
        // Calcolo errore infinito: norm(xnew - xold, inf)
        double diffNorm = 0;
        for(int i=0; i<n; i++) {
            diffNorm = Math.max(diffNorm, Math.abs(xNew.data[i] - xOld.data[i]));
        }
        
        if (diffNorm <= tol) break;
        if (Double.isNaN(diffNorm) || Double.isInfinite(diffNorm)) {
            System.out.println("Il metodo diverge! Errore NaN a iterazione " + nit);
            break;
        }

        xOld = xNew.copy();



// DENTRO IL WHILE:
temp = b.copy(); // Copia i valori di b in temp senza creare nuovi oggetti
DMatrixRMaj Bx = new DMatrixRMaj(b.numRows, 1); // Inizializza fuori dal while
// DENTRO IL WHILE:
CommonOps_DSCC.mult(B, xOld, Bx);
CommonOps_DDRM.subtract(b, Bx, temp);
// Ora temp contiene esattamente: b + (-1.0 * B * xOld)

        // RISOLUZIONE TRIANGOLARE OTTIMIZZATA (Column-based)
        // Invece di chiamare solveTriangInf, lo facciamo qui per usare 'diags' pre-calcolato
        double[] tData = temp.data;
        for (int j = 0; j < n; j++) {
            tData[j] /= diags[j];
            
            int start = L.col_idx[j];
            int end = L.col_idx[j+1];
            for (int k = start; k < end; k++) {
                int row = L.nz_rows[k];
                if (row > j) {
                    tData[row] -= L.nz_values[k] * tData[j];
                }
            }
        }
        xNew = temp.copy();
        nit++;
    }

    long durata = System.nanoTime() - inizio;
    double erroreFinale = NormOps_DDRM.normPInf(temp); // Semplificazione per l'errore

    return new Risultato(xNew, nit, durata / 1e9, erroreFinale);
}
}