/**
 * @file main.cpp
 * @brief Entry point with enhanced scientific logging
 */

#include "config.h"
#include "data_loader.h"
#include "spatial_index.h"
#include "neighborhood_mgr.h"
#include "miner.h"
#include "types.h"
#include "utils.h"
#include <iostream>
#include <fstream>
#include <chrono>
#include <windows.h>
#include <iomanip>
#include <psapi.h>

 // Helper for separating sections
void printSectionHeader(const std::string& title) {
    std::cout << "\n" << std::string(60, '=') << "\n";
    std::cout << " " << title << "\n";
    std::cout << std::string(60, '=') << "\n";
}

int main(int argc, char* argv[]) {
    auto programStart = std::chrono::high_resolution_clock::now();

    // ========================================================================
    // Step 1: Load Configuration
    // ========================================================================
    printSectionHeader("STEP 1: CONFIGURATION");
    std::string config_path = (argc > 1) ? argv[1] : "./config/config.txt";
    AppConfig config = ConfigLoader::load(config_path);

    std::cout << std::left << std::setw(25) << "Dataset Path:" << config.datasetPath << "\n";
    std::cout << std::left << std::setw(25) << "Neighbor Distance (d):" << config.neighborDistance << "\n";
    std::cout << std::left << std::setw(25) << "Min Prevalence:" << config.minPrev << "\n";

    // ========================================================================
    // Step 2: Load Data
    // ========================================================================
    printSectionHeader("STEP 2: LOAD DATA");
    const auto loadStartTime = std::chrono::high_resolution_clock::now();
    const auto instances = DataLoader::load_csv(config.datasetPath);
    const auto loadEndTime = std::chrono::high_resolution_clock::now();
    const double load_time = std::chrono::duration<double, std::milli>(loadEndTime - loadStartTime).count();

    std::cout << "[DATA] Loaded " << instances.size() << " instances.\n";
    std::cout << "[TIME] Load time: " << std::fixed << std::setprecision(2) << load_time << " ms\n";

    // ========================================================================
    // Step 3: Build Spatial Index
    // ========================================================================
    printSectionHeader("STEP 3: SPATIAL INDEXING");
    const auto indexStartTime = std::chrono::high_resolution_clock::now();
    SpatialIndex spatial_idx(config.neighborDistance);
    const auto neighborPairs = spatial_idx.findNeighborPair(instances);
    const auto indexEndTime = std::chrono::high_resolution_clock::now();
    const double idx_time = std::chrono::duration<double, std::milli>(indexEndTime - indexStartTime).count();

    std::cout << "[DATA] Found " << neighborPairs.size() << " neighbor pairs.\n";
    std::cout << "[TIME] Indexing time: " << std::fixed << std::setprecision(2) << idx_time << " ms\n";

    // ========================================================================
    // Step 4: Materialize Neighborhoods
    // ========================================================================
    printSectionHeader("STEP 4: NEIGHBORHOOD MATERIALIZATION");
    const auto materializationStartTime = std::chrono::high_resolution_clock::now();
    const std::map<FeatureType, int> featureCount = countInstancesByFeature(instances);

    NeighborhoodMgr neighbor_mgr;
    neighbor_mgr.buildFromPairs(neighborPairs, featureCount);

    NRTree orderedNRTree;
    orderedNRTree.build(neighbor_mgr, featureCount);

    const auto materializationEndTime = std::chrono::high_resolution_clock::now();
    const double mat_time = std::chrono::duration<double, std::milli>(materializationEndTime - materializationStartTime).count();

    std::cout << "[INFO] Feature Counts:\n";
    for (const auto& pair : featureCount) {
        std::cout << "   - " << std::left << std::setw(5) << pair.first << ": " << pair.second << " instances\n";
    }
    std::cout << "[TIME] Materialization time: " << std::fixed << std::setprecision(2) << mat_time << " ms\n";

    // ========================================================================
    // Step 5: Mine Colocation Patterns
    // ========================================================================
    printSectionHeader("STEP 5: MINING PROCESS");
    JoinlessMiner miner;

    // Callback đơn giản hơn, không dùng \r để tránh mất log debug
    auto progressCallback = [](int currentStep, int totalSteps, const std::string& message, double percentage) {
        // Chỉ in các mốc quan trọng hoặc message cụ thể nếu cần, 
        // ở đây ta để hàm mineColocations tự in log chi tiết nên callback có thể để trống hoặc in tối giản
        // std::cout << "[PROG] " << message << "\n"; 
        };

    auto colocations = miner.mineColocations(config.minPrev, orderedNRTree, instances, featureCount, progressCallback);

    // ========================================================================
    // Final Report
    // ========================================================================
    printSectionHeader("FINAL SUMMARY");

    const auto programEnd = std::chrono::high_resolution_clock::now();
    const double totalTimeSec = std::chrono::duration<double>(programEnd - programStart).count();
    const double maxMemory = getMemoryUsageMB();

    std::cout << std::left << std::setw(35) << "Total Prevalent Patterns:" << colocations.size() << "\n";
    std::cout << std::left << std::setw(35) << "Total Execution Time:" << std::fixed << std::setprecision(4) << totalTimeSec << " s\n";
    std::cout << std::left << std::setw(35) << "Peak Memory Usage:" << std::fixed << std::setprecision(2) << maxMemory << " MB\n";

    if (!colocations.empty()) {
        std::cout << "\n[PATTERNS FOUND]\n";
        int idx = 1;
        for (const auto& col : colocations) {
            std::cout << std::right << std::setw(3) << idx++ << ". {";
            for (size_t i = 0; i < col.size(); ++i) {
                std::cout << (i > 0 ? ", " : "") << col[i];
            }
            std::cout << "}\n";
        }
    }
    else {
        std::cout << "\n[RESULT] No patterns found satisfying the threshold.\n";
    }

    std::cout << "\nMining completed.\n";
    return 0;
}