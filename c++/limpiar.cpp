#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>

// Función para leer y procesar un archivo FASTA
std::string concatenateAndCleanFASTA(const std::string& filePath) {
    std::ifstream file(filePath);
    std::string line, concatenatedSequence;

    while (std::getline(file, line)) {
        if (line[0] != '>') {
            // Si la línea no es un encabezado, la concatenamos
            concatenatedSequence += line;
        }
    }

    // Eliminamos los caracteres 'N' de la secuencia concatenada
    concatenatedSequence.erase(std::remove(concatenatedSequence.begin(), concatenatedSequence.end(), 'N'), concatenatedSequence.end());
    return concatenatedSequence;
}

int main() {
    std::string filePath = "G5.fna"; // Archivo de entrada
    std::string cleanedSequence = concatenateAndCleanFASTA(filePath);

    // Guardamos la secuencia concatenada en un archivo de salida
    std::ofstream outputFile("G5L.fna");
    outputFile << ">Secuencia de genomas concatenada\n";
    outputFile << cleanedSequence << "\n";

    std::cout << "Se guardo correctamente'" << std::endl;
    std::cout << "Largo de la secuencia: " << cleanedSequence.size() << " caracteres" << std::endl;

    return 0;
}
