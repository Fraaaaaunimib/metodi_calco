import org.ejml.data.DMatrixSparseCSC;
import java.io.File;
import java.util.Scanner;
import java.io.IOException;

public class Main {
    public static void main(String[] args) {
        caritaMatriceManualmente("spa1.mtx");
    }

    public static void caricaMatriceManualmente(String percorso) throws IOException {
        File file = new File(percorso);
        Scanner scanner = new Scanner(file);
        String stringa = scanner.nextLine();
        while(stringa.startswith('%') || stringa.isEmpty()) {
            stringa = scanner.nextLine();
        }
        int numRighe = scanner.nextInt();
        int numColonne = scanner.nextInt();
        int nnz = scanner.nextInt();
        System.out.println("Righe: " + numRighe);
        System.out.println("Colonne: " + numColonne);
        System.out.println("Non zero: " + nnz);
        
    } 
}