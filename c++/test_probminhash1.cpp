#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>
#include <vector>
#include <cmath>
#include <random>
#include <algorithm>
#include <limits>
#include <functional>

// Se incluye minhash.hpp y las demás dependencias del repositorio del paper
#include "minhash.hpp"
#include "bitstream_random.hpp"
#include "exponential_distribution.hpp"

static const int K = 30;   // longitud del k-mer
static const uint32_t M = 30; // tamaño de la firma

// Lee los archivo FASTA y concatena las secuencias (omitiendo cabeceras)
std::string readGenome(const std::string &filename) {
    std::ifstream in(filename);
    if(!in) {
        std::cerr << "No se pudo abrir " << filename << "\n";
        exit(1);
    }
    std::string line, genome;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '>') continue;
        genome += line;
    }
    return genome;
}

// Extrae k-mers con sus frecuencias (que consideraremos como pesos)
std::unordered_map<std::string, double> extractWeightedKmers(const std::string &genome, int k=K) {
    std::unordered_map<std::string, double> kmers;
    size_t len = genome.size();
    if (len < (size_t)k) return kmers;
    for (size_t i = 0; i + k <= len; i++) {
        std::string kmer = genome.substr(i, k);
        kmers[kmer] += 1.0;
    }
    return kmers;
}

// Definimos el tipo de elemento que pasaremos al minhash:
// Se usa un vector de pares (D, double), donde D es el k-mer hasheado (uint64_t) y double es su peso.
struct KmerItem {
    uint64_t element; // hash del k-mer
    double weight;
};

// Dado un KmerItem, extrae el identificador D (aquí uint64_t)
struct ExtractFunction {
    uint64_t operator()(const KmerItem &item) const {
        return item.element;
    }
};

// weightFunction: Dado un KmerItem, retorna el peso
struct WeightFunction {
    double operator()(const KmerItem &item) const {
        return item.weight;
    }
};

// Dado un element (uint64_t), genera un WyrandBitStream (pára la generación de numeros aleatorios)
struct RngFunction {
    WyrandBitStream operator()(uint64_t element) const {
        // Añadir una semilla global fija diferente o combinar el elemento con 
        // un seed fijo distinto para asegurar que cada elemento genere secuencias 
        // aleatorias reproducibles pero únicas.
        return WyrandBitStream(element, 123456789ULL);
    }
};


// Calcula la similitud Jaccard ponderada aproximada: cuenta cuántos componentes son iguales.
double estimateJaccard(const std::vector<uint64_t> &sigA, const std::vector<uint64_t> &sigB) {
    if (sigA.size() != sigB.size()) {
        std::cerr << "Las firmas tienen tamaños diferentes\n";
        return 0.0;
    }
    int m = (int)sigA.size();
    int count = 0;
    for (int i=0; i<m; i++) {
        if (sigA[i] == sigB[i]) count++;
    }
    return double(count)/m;
}

int main() {
    std::vector<std::string> files = {"G1L.fna","G2L.fna","G3L.fna","G4L.fna","G5L.fna"};

    // Hash para k-mers
    std::hash<std::string> str_hash;

    // Cargamos y extraemos k-mers ponderados
    std::vector<std::vector<KmerItem>> allWeightedSets;
    for (auto &f : files) {
        std::string genome = readGenome(f);
        auto wkmers = extractWeightedKmers(genome, K);

        // Convertimos a vector de KmerItem
        std::vector<KmerItem> items;
        items.reserve(wkmers.size());
        for (auto &kv : wkmers) {
            uint64_t h = str_hash(kv.first);
            KmerItem it {h, kv.second};
            items.push_back(it);
        }
        allWeightedSets.push_back(std::move(items));
    }

    // Ejecutamos el algoritmo ProbMinHash1 con parámetros
    // D = uint64_t
    // E = ExtractFunction
    // R = RngFunction
    // W = WeightFunction
    ProbMinHash1<uint64_t, ExtractFunction, RngFunction, WeightFunction> pmh(M, ExtractFunction(), RngFunction(), WeightFunction());

    // Generamos las firmas
    std::vector<std::vector<uint64_t>> signatures;
    for (auto &dataset : allWeightedSets) {
        auto result = pmh(dataset); 
        // result es un vector<D> con los elementos que definieron cada componente
        // Lo convertimos a vector<uint64_t> representando la firma.
        signatures.push_back(result);
    }

    // Calculamos la similitud entre todos los pares
    int n = (int)signatures.size();
    for (int i=0; i<n; i++) {
        for (int j=i+1; j<n; j++) {
            double sim = estimateJaccard(signatures[i], signatures[j]);
            std::cout << "Jaccard ponderada entre " << files[i] << " y " << files[j] << ": " << sim << "\n";
        }
    }

    return 0;
}
