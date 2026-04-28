package com.funzioni;
import org.ejml.data.DMatrixRMaj;
import org.ejml.data.DMatrixSparseCSC;
import org.ejml.sparse.csc.CommonOps_DSCC;

import java.io.File;
import java.util.Scanner;
import java.io.IOException;


public class Funzioni {
    public static DMatrixSparseCSC caricaMatriceManualmente(String percorso) { // Metodo per caricare una matrice da un file .mtx
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
        System.out.println("Righe: " + numRighe);
        System.out.println("Colonne: " + numColonne);
        System.out.println("Non zero: " + nnz);

        DMatrixSparseCSC matrice = new DMatrixSparseCSC(numRighe, numColonne, nnz); // Crea una matrice sparsa con il numero di righe, colonne e non zero specificati
        for(int i = 0; i < nnz; i++) {
            int righe = scanner.nextInt() -1; // Sottrai 1 per convertire da indice 1-based a 0-based
            int colonne = scanner.nextInt() -1; // Sottrai 1 per convertire da indice 1-based a 0-based
            double valore = scanner.nextDouble();
            matrice.set(righe, colonne, valore);
        }

        scanner.close();
        
        return matrice;
    }

    public static boolean diagonaleZeri(DMatrixSparseCSC A) {
        for(int i = 0; i < A.numRows; i++) {
            if(A.get(i, i) == 0.0){
                return false;
            }
        }
        return true;
    }


    public static double normaInfinito(double[] v, double[] w){
        double max = 0.0;
        for(int i = 0; i < v.length; i++) {
            double diff = Math.abs(v[i] - w[i]);
            if(diff > max) {
                max = diff;
            }
        }
        return max;
    }

    public static DMatrixSparseCSC triang_inf(final DMatrixSparseCSC A){
    DMatrixSparseCSC aT = A.copy();
    for(int i = 0; i < A.numRows; i++){
        for(int j = 0; j < A.numCols; j++){
            if(j > i){ // Se l'indice di colonna è maggiore di quello di riga
                aT.set(i, j, 0.0); // Azzera l'elemento sopra la diagonale
            }
        }
    }
    return aT;
    }

    public static void solveTriangInf(DMatrixSparseCSC L, DMatrixRMaj B) {
    int n = L.numCols;
    double[] bData = B.data;

    // Sostituzione in avanti ottimizzata per il formato a colonne (CSC)
    for (int j = 0; j < n; j++) {
        // 1. Troviamo l'elemento sulla diagonale L(j,j)
        double diag = 0;
        int start = L.col_idx[j];
        int end = L.col_idx[j+1];

        // Cerchiamo la diagonale nella colonna j
        for (int k = start; k < end; k++) {
            if (L.nz_rows[k] == j) {
                diag = L.nz_values[k];
                break;
            }
        }

        if (Math.abs(diag) < 1e-15) continue;

        // 2. Calcoliamo il valore definitivo dell'incognita j
        bData[j] /= diag;

        // 3. Aggiorniamo (sottraiamo) l'effetto di x(j) su tutte le righe sottostanti
        // Questo è il segreto: scorriamo la COLONNA, che in CSC è velocissimo!
        for (int k = start; k < end; k++) {
            int riga = L.nz_rows[k];
            if (riga > j) {
                bData[riga] -= L.nz_values[k] * bData[j];
            }
        }
    }
}
public static void ritornoValori(double tol, Risultato r){
            System.out.print("| Tol = " + tol);
            System.out.print(" | Iterazioni = " + r.getNit());
            System.out.print(" | Errore = " + r.getErrore());
            System.out.println(" | Tempo = " + r.getTempo() + " | ");
            System.out.println("----------------------");
    }
}

