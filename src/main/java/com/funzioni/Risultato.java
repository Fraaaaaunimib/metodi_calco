package com.funzioni;

import org.ejml.data.DMatrixRMaj;

public class Risultato {
    private DMatrixRMaj x;
    private int nit;
    private double tempo;
    private double errore;

    public Risultato(DMatrixRMaj x, int nit, double tempo, double errore) {
        this.x = x;
        this.nit = nit;
        this.tempo = tempo;
        this.errore = errore;
    }

    


    public DMatrixRMaj getX() {
        return x;
    }

    public void setX(DMatrixRMaj x) {
        this.x = x;
    }

    public int getNit() {
        return nit;
    }

    public void setNit(int nit) {
        this.nit = nit;
    }

    public double getTempo() {
        return tempo;
    }

    public void setTempo(double tempo) {
        this.tempo = tempo;
    }

    public double getErrore() {
        return errore;
    }

    public void setErrore(double errore) {
        this.errore = errore;
    }

}
