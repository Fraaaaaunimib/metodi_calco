package com.funzioni;

// Risultato dell'esecuzione di un metodo iterativo
public class Risultato {
    private int nit; // numero di iterazioni
    private double tempo; // tempo di esecuzione
    private double errore; // errore relativo
    
    /*
    Costruttore pubblico da usare alla fine di ciascun metodo iterativo
    Parametri: nit (numero di iterazioni), tempo (tempo di esecuzione), errore (errore relativo)
    */
    public Risultato(int nit, double tempo, double errore) {
        this.nit = nit;
        this.tempo = tempo;
        this.errore = errore;
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
