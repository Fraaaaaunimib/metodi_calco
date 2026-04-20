import org.ejml.data.DMatrixSparseCSC;
import org.ejml.ops.MatrixIO;
import java.io.IOException;

public class Main {
    public static void main(String[] args) {
        try {
            // Legge il file .mtx e lo carica in una matrice sparsa
            // DMatrixSparseCSC è il formato ideale per matrici con molti zeri
            DMatrixSparseCSC matrice = MatrixIO.loadMatrixMarketDCSC("/workspaces/metodi_calco/spa1.mtx");

            System.out.println("Matrice caricata con successo!");
            System.out.println("Righe: " + matrice.numRows);
            System.out.println("Colonne: " + matrice.numCols);

        } catch (IOException e) {
            System.err.println("Errore durante la lettura del file: " + e.getMessage());
        }
    }
}