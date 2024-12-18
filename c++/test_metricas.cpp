#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <chrono>
#include <cstdlib>
#include <cstring>
#include <unordered_map>

// Función para medir tiempo de ejecución
double measureExecutionTime(const std::string &command) {
    auto start = std::chrono::high_resolution_clock::now();
    int ret = std::system(command.c_str());
    (void)ret; // ignorar valor de retorno
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    return elapsed.count();
}

// Función para medir uso de memoria 
size_t getMemoryUsage() {
    std::ifstream statusFile("/proc/self/status");
    std::string line;
    while(std::getline(statusFile, line)) {
        if(line.rfind("VmRSS:", 0) == 0) {
            std::string val;
            for (int i = 6; i < (int)line.size(); i++) {
                if (std::isdigit((unsigned char)line[i])) {
                    val.push_back(line[i]);
                } else if (!val.empty()) {
                    break;
                }
            }
            if(!val.empty()) {
                return std::stoul(val) * 1024; // transformar kB en bytes
            }
        }
    }
    return 0;
}

// Función que ejecuta el comando y mide memoria antes y después.
// Devuelve un par <tiempo, memoria>
std::pair<double, size_t> measureTimeAndMemory(const std::string &command) {
    size_t mem_before = getMemoryUsage();
    double time = measureExecutionTime(command);
    size_t mem_after = getMemoryUsage();
    size_t mem_used = (mem_after > mem_before) ? (mem_after - mem_before) : 0;
    return {time, mem_used};
}

int main() {
    std::string minhash_cmd = "./test_minhash";
    std::string probminhash_cmd = "./test_probminhash1";

    std::cout << "=== MinHash Metrics ===\n";
    auto [time_min, mem_min] = measureTimeAndMemory(minhash_cmd);
    std::cout << "Tiempo de ejecución: " << time_min << " s\n";
    std::cout << "Uso de memoria: " << mem_min << " bytes\n";

    std::cout << "\n=== ProbMinHash Metrics ===\n";
    auto [time_prob, mem_prob] = measureTimeAndMemory(probminhash_cmd);
    std::cout << "Tiempo de ejecución: " << time_prob << " s\n";
    std::cout << "Uso de memoria: " << mem_prob << " bytes\n";

    return 0;
}
