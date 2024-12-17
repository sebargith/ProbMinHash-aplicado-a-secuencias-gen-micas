#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <functional>
#include <limits>
#include <random>
#include <algorithm>

static const int K = 30;     // longitud del k-mer
static const uint32_t M = 128; // tamaño de la firma MinHash 

// Lee un archivo FASTA concatenando todas las líneas (omitiendo cabeceras)
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

// Extrae k-mers y retorna un set de k-mers únicos
std::unordered_set<std::string> extractUniqueKmers(const std::string &genome, int k=K) {
    std::unordered_set<std::string> kmers;
    size_t len = genome.size();
    if (len < (size_t)k) return kmers;
    for (size_t i=0; i + k <= len; i++) {
        std::string kmer = genome.substr(i, k);
        kmers.insert(kmer);
    }
    return kmers;
}

// Crea M "semillas" diferentes para generar M funciones hash.
std::vector<uint64_t> generateHashSeeds(uint32_t m) {
    std::mt19937_64 gen(12345); // semilla fija para reproducibilidad
    std::uniform_int_distribution<uint64_t> dist;
    std::vector<uint64_t> seeds(m);
    for (uint32_t i=0; i<m; i++) {
        seeds[i] = dist(gen);
    }
    return seeds;
}

// Calculamos la firma MinHash de un conjunto de k-mers dado
// Usamos std::hash<std::string> combinado con una semilla diferente para cada componente.
std::vector<uint64_t> computeMinHashSignature(const std::unordered_set<std::string> &kmers, 
                                              const std::vector<uint64_t> &seeds) {
    std::hash<std::string> str_hash;
    uint32_t m = (uint32_t)seeds.size();
    std::vector<uint64_t> signature(m, std::numeric_limits<uint64_t>::max());

    for (auto &kmer : kmers) {
        uint64_t baseHash = str_hash(kmer);
        // combinamos con la seed y tomamos el mínimo.
        for (uint32_t i=0; i<m; i++) {
            uint64_t h = baseHash ^ seeds[i];
            if (h < signature[i]) {
                signature[i] = h;
            }
        }
    }

    return signature;
}

// Calculamos la similitud Jaccard estimada a partir de dos firmas minhash
double estimateJaccardFromMinHash(const std::vector<uint64_t> &sigA, const std::vector<uint64_t> &sigB) {
    if (sigA.size() != sigB.size()) {
        std::cerr << "Firmas con tamaños distintos, no se puede calcular similitud\n";
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
    // Archivos de genomas
    std::vector<std::string> files = {"G1L.fna","G2L.fna","G3L.fna","G4L.fna","G5L.fna"};

    // Leer genomas y extraer k-mers
    std::vector<std::unordered_set<std::string>> allSets;
    for (auto &f : files) {
        std::string genome = readGenome(f);
        auto kmers = extractUniqueKmers(genome, K);
        allSets.push_back(std::move(kmers));
    }

    // Generamos las seeds para minhash
    auto seeds = generateHashSeeds(M);

    // Calculamos las firmas minhash
    std::vector<std::vector<uint64_t>> signatures;
    for (auto &kmerset : allSets) {
        auto sig = computeMinHashSignature(kmerset, seeds);
        signatures.push_back(sig);
    }

    // Calculamos similitud entre todos los pares
    int n = (int)signatures.size();
    for (int i=0; i<n; i++) {
        for (int j=i+1; j<n; j++) {
            double sim = estimateJaccardFromMinHash(signatures[i], signatures[j]);
            std::cout << "Similitud Jaccard (MinHash) entre " << files[i] << " y " << files[j] << ": " << sim << "\n";
        }
    }

    return 0;
}
