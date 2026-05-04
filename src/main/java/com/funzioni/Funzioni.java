package com.funzioni;
import org.ejml.data.DMatrixRMaj;
import org.ejml.data.DMatrixSparseCSC;
import org.ejml.dense.row.CommonOps_DDRM;
import org.ejml.dense.row.NormOps_DDRM;
import org.ejml.sparse.csc.CommonOps_DSCC;

import java.io.File;
import java.util.Scanner;
import java.io.IOException;

// Metodi di supporto ai metodi iterativi
public class Funzioni {
    
    /*
    Caricamento matrice da File .mtx
    Parametri: percorso (percorso del file .mtx)
    Ritorno: matrice sparsa
    */
    public static DMatrixSparseCSC caricaMatriceManualmente(String percorso) {
        File file = new File(percorso);
        Scanner scanner;
        try{
            scanner = new Scanner(file);
        }
        catch(IOException e) {
            System.out.println("Errore nella lettura del file: " + e.getMessage());
            return null;
        }
        String stringa = scanner.next();
        while(stringa.startsWith("%") || stringa.isEmpty()) {
            stringa = scanner.nextLine();
        }
        int numRighe = scanner.nextInt();
        int numColonne = scanner.nextInt();
        int nnz = scanner.nextInt();
        
        DMatrixSparseCSC matrice = new DMatrixSparseCSC(numRighe, numColonne, nnz);
        for(int i = 0; i < nnz; i++) {
            int righe = scanner.nextInt() -1; 
            int colonne = scanner.nextInt() -1; 
            double valore = scanner.nextDouble();
            matrice.set(righe, colonne, valore);
        }
        
        scanner.close();
        
        return matrice;
    }
    
    /*
    Controlla esistenza di zeri sulla diagonale
    Parametri: A (matrice sparsa)
    Ritorno: true se non ha zeri sulla diagonale, false altrimenti 
    */
    public static boolean isNotDiagonaleZero(DMatrixSparseCSC A) {
        for(int i = 0; i < A.numRows; i++) {
            if(A.get(i, i) == 0.0){
                return false;
            }
        }
        return true;
    }
    
    /*
    Crea la matrice triangolare inferiore a partire da una matrice sparsa
    Parametri: A (matrice sparsa)
    Ritorno: matrice triangolare creata a partire dal parametro
    */
    public static DMatrixSparseCSC triang_inf(final DMatrixSparseCSC A) {
        DMatrixSparseCSC L = new DMatrixSparseCSC(A.numRows, A.numCols, A.nz_length);
        for (int col = 0; col < A.numCols; col++) {
            int idx0 = A.col_idx[col];
            int idx1 = A.col_idx[col + 1];
            for (int idx = idx0; idx < idx1; idx++) {
                int row = A.nz_rows[idx];
                if (row >= col) {
                    L.set(row, col, A.nz_values[idx]);
                }
            }
        }
        return L;
    }
    
    /*
    Risolve il sistema Lx=b usando la forward substitution
    Parametri: L (matrice sparsa triangolare inferiore) e b (vettore dei termini noti)
    Ritorno: x (vettore soluzione)
    */
    public static DMatrixRMaj solveTriangInf(DMatrixSparseCSC L, DMatrixRMaj b) {
        
        int n = L.numRows;
        
        if (L.numCols != n) {
            throw new IllegalArgumentException("Matrice non quadrata");
        }
        
        DMatrixRMaj x = new DMatrixRMaj(n, 1);
        
        for (int i = 0; i < n; i++) {
            
            double somma = 0.0;
            
            for (int j = 0; j < i; j++) {
                somma += L.get(i, j) * x.get(j, 0);
            }
            
            double diag = L.get(i, i);
            
            if (diag == 0.0) {
                throw new ArithmeticException("Zero sulla diagonale in riga " + i);
            }
            
            double valore = (b.get(i, 0) - somma) / diag;
            x.set(i, 0, valore);
        }
        
        return x;
    }
    
    /*
    Stampa le soluzioni di un'iterazione di uno dei metodi iterativi
    Parametri: tol (tolleranza) e r (Risultato dell'iterazione)
    */
    public static void ritornoValori(double tol, Risultato r){
        System.out.print("| Tol = " + tol);
        System.out.print(" | Iterazioni = " + r.getNit());
        System.out.print(" | Errore = " + r.getErrore());
        System.out.println(" | Tempo = " + r.getTempo() + " | ");
        System.out.println("----------------------");
    }
    
    /*
    Calcolo di ||A * x0 - b||_2 / ||b||_2
    Parametri: A (matrice sparsa), x0 (vettore della soluzione del metodo), b (vettore dei termini noti)
    Ritorno: risultato della frazione
    */
    public static double calcoloNormaDiff(final DMatrixSparseCSC A, final DMatrixRMaj x0, final DMatrixRMaj b) {
        DMatrixRMaj workspace = new DMatrixRMaj(A.numCols, 1);
        
        CommonOps_DSCC.mult(A, x0, workspace);
        CommonOps_DDRM.subtract(workspace, b, workspace);
        double normaDiff = NormOps_DDRM.normP2(workspace);
        double normaInf = NormOps_DDRM.normP2(b);
        normaDiff = normaDiff / normaInf;
        
        return normaDiff;
    }

    /*
    Calcolo dell'errore relativo valido per tutti i metodi
    Parametri: n (numero di righe della matrice iniziale), xNew (vettore dei risultati del metodo)
    Ritorno: errore relativo
    */
    public static double calcoloErroreRelativo(int n, DMatrixRMaj xNew){
        DMatrixRMaj xEsatta = new DMatrixRMaj(n, 1);
        xEsatta.fill(1.0);

        // Calcola differenza xNew - xEsatta
        DMatrixRMaj diffErr = new DMatrixRMaj(n, 1);
        CommonOps_DDRM.subtract(xNew, xEsatta, diffErr);

        // Errore relativo
        return NormOps_DDRM.normP2(diffErr) / NormOps_DDRM.normP2(xEsatta);
    }

    /*
    Controlla che la matrice sia quadrata e che la dimensione del vettore x0 sia uguale alla matrice
    Parametri: n (numero di colonne della matrice), m (numero di righe della matrice), L (numero di righe del vettore x0)
    Nota: lancia eccezioni per i controlli
    */
    public static void controlloDimensione(int n, int m, int L) throws Exception{
        if (n != m) {
            throw new Exception("La matrice A deve essere quadrata.");
        } else if (L != m) {
            throw new Exception("La dimensione di A deve essere uguale a x0.");
        } else {
            return;
        }
    }
}

