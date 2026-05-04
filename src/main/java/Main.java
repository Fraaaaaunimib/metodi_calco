import org.ejml.data.DMatrixRMaj;
import org.ejml.data.DMatrixSparseCSC;
import org.ejml.sparse.csc.CommonOps_DSCC;

import com.funzioni.*;
import com.metodi.Gauss;
import com.metodi.Jacobi;
import com.metodi.gradient;
import com.metodi.gradientC;

public class Main {

    static int maxIter = 20000;
    
    public static void main(String[] args ) {

        if(args.length == 0){
            System.err.println("Non hai messo nessun argoment");
            return;
        }
        else if(args.length == 1){
            System.err.println("Hai messo un solo argomento, ne servono 2");
            return;
        }

        String nomeMatrice = args[0];
        double tol = 0.0;
        try{
            tol = Double.parseDouble(args[1]);
        } catch(NumberFormatException e){
            System.err.println("Eccezione: " + e);
            return;
        }
        DMatrixSparseCSC matrice = Funzioni.caricaMatriceManualmente(nomeMatrice);
        
        DMatrixRMaj x0 = new DMatrixRMaj(matrice.numRows, 1);
            for(int i = 0; i < x0.getNumRows(); i++) {
                x0.set(i, 0.0);
            }

            DMatrixRMaj x = new DMatrixRMaj(matrice.numRows, 1);
            for(int i = 0; i < x.getNumRows(); i++) {
                x.set(i, 1.0);
            }

            DMatrixRMaj b = new DMatrixRMaj(matrice.numRows, 1);
            CommonOps_DSCC.mult(matrice, x, b);

        System.out.println("Metodo di Jacobi");
        Risultato rJacobi = Jacobi.jacobi(matrice, b, x0, maxIter, tol);
        Funzioni.ritornoValori(tol, rJacobi);
        System.out.println();
        System.out.println("Metodo di Gauss");
        Risultato rGauss = Gauss.gauss(matrice, b, x0, maxIter, tol);
        Funzioni.ritornoValori(tol, rGauss);
        System.out.println();
        System.out.println("Metodo del Gradiente");
        Risultato rGradiente = gradient.gradiente(matrice, b, x0, maxIter, tol);
        Funzioni.ritornoValori(tol, rGradiente);
        System.out.println();
        System.out.println("Metodo del Gradiente Coniugato");
        Risultato rGradienteC = gradientC.gradienteC(matrice, b, x0, maxIter, tol);
        Funzioni.ritornoValori(tol, rGradienteC);
    }
}