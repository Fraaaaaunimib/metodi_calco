import org.ejml.data.DMatrixSparseCSC;
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
        
        System.out.println("Caricamento matrice:");

        DMatrixSparseCSC matrice = new DMatrixSparseCSC(numRighe, numColonne, nnz); // Crea una matrice sparsa con il numero di righe, colonne e non zero specificati
        for(int i = 0; i < nnz; i++) {
            int righe = scanner.nextInt() -1; // Sottrai 1 per convertire da indice 1-based a 0-based
            int colonne = scanner.nextInt() -1; // Sottrai 1 per convertire da indice 1-based a 0-based
            double valore = scanner.nextDouble();
            matrice.set(righe, colonne, valore);
        }
        return matrice;
        /*for(int i = 0; i < 5; i++) { // Stampa i primi 5x4 elementi non zero della matrice
            for(int j = 0; j < 4; j++) {
                double valore = matrice.get(i, j);
                if(valore != 0) {
                    System.out.println("Elemento (" + i + ", " + j + "): " + valore);
                }
            }
        }*/
    } 

    public static boolean diagonaleZeri(DMatrixSparseCSC A) {
        for(int i = 0; i < A.numRows; i++) {
            if(A.get(i, i) == 0.0){
                return false;
            }
        }
        return true;
    }

    

}
