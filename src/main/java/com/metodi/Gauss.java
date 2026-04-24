package com.metodi;
import com.funzioni.*;

import java.util.Arrays;

import org.ejml.data.DMatrix;
import org.ejml.data.DMatrixRMaj;
import org.ejml.data.DMatrixSparseCSC;
import org.ejml.dense.row.CommonOps_DDRM;
import org.ejml.dense.row.NormOps_DDRM;
import org.ejml.sparse.csc.CommonOps_DSCC;

public class Gauss{
    public Risultato gauss(DMatrixSparseCSC A, DMatrixRMaj b, DMatrixRMaj x0, int maxIter, double tol) {
        int m = A.numRows;
        int n = A.numCols;
        int L = x0.getNumRows();

        if(m!= n){
            System.out.println("La matrice A deve essere quadrata.");
            return null;
        }
        else if(L!=m){
            System.out.println("La dimensione di A deve essere uguale alla dimensione di x0.");
            return null;
        }

        DMatrixSparseCSC A_triang = Funzioni.triang_inf(A);
        DMatrixSparseCSC B = new DMatrixSparseCSC(A.numRows, A.numCols);
        CommonOps_DSCC.add(1.0, A, -1.0, A_triang, B, null, null); // Calcola B = A - A_triang

        DMatrixRMaj xOld = x0.copy();
        DMatrixRMaj xNew = xOld.copy();

        DMatrixRMaj diffOld = new DMatrixRMaj(xOld.getNumRows(), 1);
        Arrays.fill(diffOld.data, Double.MAX_VALUE);

       // CommonOps_DDRM.subtract(DMatrixRMaj.wrap(xNew.getNumRows(), 1, xNew.getData()), xOld, diffOld);

        int nit = 0;

        long inizio = System.nanoTime();
        DMatrixRMaj temp = new DMatrixRMaj(A.numRows, 1);


        while(NormOps_DDRM.normPInf(diffOld) > tol && nit < maxIter){
            xOld = xNew.copy();
            CommonOps_DSCC.mult(B, xOld, temp);
            CommonOps_DDRM.subtract(b,temp, temp);

        }
        
    }
}