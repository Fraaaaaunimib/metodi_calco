# Progetto 1-Alternativo
Metodi del Calcolo Scientifico 2025/2026

Cattaneo Francesco (matricola 900411), Carpio Herreros Marco (matricola 899802)

## Build
Per la prima esecuzione: `mvn clean compile` nella root della repository.

## Esecuzione
Il programma si aspetta due argomenti: `args[0]` (nome del file della matrice .mtx) e `args[1]` (tolleranza).

Il comando è: `mvn exec:java -Dexec.args="percorso_relativo/nome_matrice.mtx 1e-6"`.
