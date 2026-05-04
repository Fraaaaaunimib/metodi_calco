package com.metodi;

import com.funzioni.*;

import org.ejml.data.DMatrixRMaj;
import org.ejml.data.DMatrixSparseCSC;
import org.ejml.dense.row.CommonOps_DDRM;
import org.ejml.sparse.csc.CommonOps_DSCC;

public class Gauss {
    
    /*
    Implementazione del metodo di Gauss-Seidel per matrici sparse
    Parametri: A (matrice sparsa di partenza), b (vettore di termini noti), x0 (vettore nullo), maxIter (numero massimo di iterazioni, deve essere >=20000), tol (tolleranza)
    Ritorno: oggetto di tipo Risultato (numero di iterazioni, tempo di esecuzione, errore relativo)
    Errore:
    - Se la matrice A non è quadrata
    - Se la lunghezza di x0 è diversa da A
    */
    public static Risultato gauss(DMatrixSparseCSC A, DMatrixRMaj b, DMatrixRMaj x0, int maxIter, double tol) {
        
        int n = A.numCols;
        int m = A.numRows;
        
        try {
            Funzioni.controlloDimensione(n, m, x0.numRows);
        } catch (Exception e){
            System.err.println(e);
        }
        
        // L = triangolare inferiore di A
        DMatrixSparseCSC L = Funzioni.triang_inf(A);
        // B = A - L
        DMatrixSparseCSC B = new DMatrixSparseCSC(n, n, A.nz_length);
        
        CommonOps_DSCC.add(1.0, A, -1.0, L, B, null, null);
        
        DMatrixRMaj xNew = x0.copy();
        DMatrixRMaj xOld = new DMatrixRMaj(n, 1);
        
        DMatrixRMaj temp = new DMatrixRMaj(n, 1);
        DMatrixRMaj diff = new DMatrixRMaj(n, 1);
        
        int nit = 0;
        
        double normaDiff = Funzioni.calcoloNormaDiff(A, xNew, b);
        
        long inizio = System.nanoTime();
        
        while (normaDiff > tol && nit < maxIter) {
            
            // xOld = xNew
            CommonOps_DDRM.scale(1.0, xNew, xOld);
            
            // diff = b - B*xOld
            CommonOps_DSCC.mult(B, xOld, temp);
            CommonOps_DDRM.subtract(b, temp, diff);
            
            // risolvi L xNew = diff
            DMatrixRMaj xNext = Funzioni.solveTriangInf(L, diff);
            
            // aggiorna soluzione
            xNew = xNext;
            
            normaDiff = Funzioni.calcoloNormaDiff(A, xNew, b);
            
            nit++;
        }
        
        long durata = System.nanoTime() - inizio;
        double tempoSecondi = durata / 1_000_000_000.0;
        
        return new Risultato(nit, tempoSecondi, Funzioni.calcoloErroreRelativo(n, xNew));
    }
}