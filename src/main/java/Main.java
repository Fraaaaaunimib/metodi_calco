import org.ejml.data.DMatrixSparseCSC;
import com.funzioni.Funzioni;
import com.metodi.Jacobi;

public class Main {
    public static void main(String[] args) {
        DMatrixSparseCSC matrice = Funzioni.caricaMatriceManualmente("spa1.mtx");

        Jacobi jacobi = new Jacobi();
        double[] b = new double[matrice.numRows];
        for(int i = 0; i < b.length; i++) {
            b[i] = Math.random()*100;
        }
        double[] x0 = new double[matrice.numRows];
        for(int i = 0; i < b.length; i++) {
            x0[i] = 0.0;
        }

        double[] jacobiR = jacobi.jacobi(matrice, b, x0, 100, 1e-6);
        System.out.println("Tempo di esecuzione: " + jacobiR[0] + " secondi");
        System.out.println("Numero di iterazioni: " + jacobiR[1]);
        System.out.println("Errore: " + jacobiR[2]);

    }

    
    
}