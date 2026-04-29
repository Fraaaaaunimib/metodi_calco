package com.funzioni;
import org.ejml.data.DMatrixRMaj;
import org.ejml.data.DMatrixSparseCSC;
import org.ejml.dense.row.CommonOps_DDRM;
import org.ejml.dense.row.NormOps_DDRM;
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

    //da rimuovere?
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

    public static DMatrixRMaj solveTriangInf(DMatrixSparseCSC L, DMatrixRMaj B) {
        double sum = 0.0; int i = 0; int j = 0;
        DMatrixRMaj X = new DMatrixRMaj(L.numRows);

        for(i = 0; i< L.numRows; i++){
            for(j = 0; j < i; j++){
                sum += L.get(i, j) * X.get(j);
            }
            X.set(i, (1.0 / L.get(i, j)) * (B.get(i) - sum));
        }

        return X;

    }
public static void ritornoValori(double tol, Risultato r){
            System.out.print("| Tol = " + tol);
            System.out.print(" | Iterazioni = " + r.getNit());
            System.out.print(" | Errore = " + r.getErrore());
            System.out.println(" | Tempo = " + r.getTempo() + " | ");
            System.out.println("----------------------");
    }

    public static double calcoloNormaDiff(final DMatrixSparseCSC A, final DMatrixRMaj x0, final DMatrixRMaj b,
            final DMatrixRMaj workspace) {
        // ||A * x0 - b||_inf / ||b||_inf
        CommonOps_DSCC.mult(A, x0, workspace);
        CommonOps_DDRM.subtract(workspace, b, workspace);
        double normaDiff = NormOps_DDRM.normP2(workspace);
        double normaInf = NormOps_DDRM.normP2(b);
        normaDiff = normaDiff / normaInf;

        return normaDiff;
    }
}

