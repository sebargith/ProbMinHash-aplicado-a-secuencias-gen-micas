# ProbMinHash para Análisis de Genomas

Sebastián Rosas Urra, Tópicos en grandes volúmenes de datos, UdeC

Este proyecto, para el curso de Tópicos en grandes volumenes de datos, implementa y compara dos algoritmos para estimar similitudes entre genomas a partir de k-mers:

- **MinHash:** Método tradicional que estima la similitud Jaccard en conjuntos no ponderados.
- **ProbMinHash1:** Extensión basada en el paper [*ProbMinHash – A Class of Locality-Sensitive Hash Algorithms for the (Probability) Jaccard Similarity* - Otmar Ertl (2022)] que soporta conjuntos ponderados y estima la similitud Jaccard probabilística.

La implementación toma como entrada secuencias genómicas (archivos `.fna`) y extrae k-mers para construir firmas que permiten estimar la similitud entre genomas.

## Características

- Extracción de k-mers a partir de secuencias genómicas.
- Cálculo de MinHash (no ponderado) y ProbMinHash1 (ponderado).
- Parámetros configurables: longitud del k-mer (K) y tamaño de la firma (M).
- Medición de tiempo de ejecución, uso de memoria y similitud estimada.
- Comparación entre ambos métodos en distintos escenarios.

## Prerrequisitos

- Entorno Linux (Ubuntu preferentemente)
- Compilador C++
- Bibliotecas estándar de C++.

## Instrucciones de Compilación

1. Clonar el repositorio
2. Compilar "test_minhash.cpp"
   g++ -o test_minhash test_minhash.cpp
3. Compilar "test_probminhash1.cpp"
   g++ -o test_probminhash1 test_probminhash1.cpp
4. (Opcional) Compilar "test_metricas.cpp
   g++ -o test_metricas test_metricas.cpp

