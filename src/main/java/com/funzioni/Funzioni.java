package com.funzioni;
import org.ejml.data.DMatrixRMaj;
import org.ejml.data.DMatrixSparseCSC;
import org.ejml.dense.row.CommonOps_DDRM;
import org.ejml.dense.row.NormOps_DDRM;
import org.ejml.interfaces.decomposition.CholeskySparseDecomposition_F64;
import org.ejml.sparse.csc.CommonOps_DSCC;
import org.ejml.sparse.csc.factory.DecompositionFactory_DSCC;
import org.ejml.sparse.csc.factory.DecompositionFactory_DSCC;
import org.ejml.interfaces.decomposition.CholeskySparseDecomposition_F64;

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

    public static DMatrixSparseCSC triang_inf(final DMatrixSparseCSC A) {
    DMatrixSparseCSC L = new DMatrixSparseCSC(A.numRows, A.numCols, A.nz_length);
    for (int col = 0; col < A.numCols; col++) {
        int idx0 = A.col_idx[col];
        int idx1 = A.col_idx[col + 1];
        for (int idx = idx0; idx < idx1; idx++) {
            int row = A.nz_rows[idx];
            if (row >= col) { // triangolare inferiore: riga >= colonna
                L.set(row, col, A.nz_values[idx]);
            }
        }
    }
    return L;
}

    /*
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
/*
    public static boolean metodoPotenze(DMatrixSparseCSC A, int nMax, double tol){
        DMatrixRMaj z = new DMatrixRMaj(A.numRows, 1);
        DMatrixRMaj q = new DMatrixRMaj(A.numRows, 1);
        q.fill(1.0);

        NormOps_DDRM.normalizeF(q);
        double lambda = 0.0;
        for(int i = 0; i < nMax; i++){
            CommonOps_DSCC.mult(A, q, z);
            double zNorma = NormOps_DDRM.normP2(z);
            DMatrixRMaj q1 = new DMatrixRMaj(A.numRows, 1);
            CommonOps_DDRM.divide(z, zNorma, q1);
            DMatrixRMaj q1T = new DMatrixRMaj(1, A.numRows);
            CommonOps_DDRM.transpose(q1T);
            CommonOps_DSCC.mult(A, q1, q1);
            lambda = CommonOps_DDRM.dot(q1, q1T);
            
        }
    }
        */
}

