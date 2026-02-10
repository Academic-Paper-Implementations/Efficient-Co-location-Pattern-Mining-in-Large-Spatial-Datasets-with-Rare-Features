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
#include <iomanip>

 //Show memmory usage
#include <windows.h>
#include <psapi.h>
#include <stdio.h>
#pragma comment(lib, "psapi.lib")

int main(int argc, char* argv[]) {
    auto programStart = std::chrono::high_resolution_clock::now();

    // ========================================================================
    // Step 1: Load Configuration
    // ========================================================================
    std::cout << "Running... (Results will be saved to result.txt)\n";
    std::string config_path = (argc > 1) ? argv[1] : "./config/config.txt";
    AppConfig config = ConfigLoader::load(config_path);

    // ========================================================================
    // Step 2: Load Data
    // ========================================================================
    auto instances = DataLoader::load_csv(config.datasetPath);

    // ========================================================================
    // Step 3: Build Spatial Index
    // ========================================================================
    SpatialIndex spatial_idx(config.neighborDistance);
    auto neighborPairs = spatial_idx.findNeighborPair(instances);

    // ========================================================================
    // Step 4: Materialize Neighborhoods
    // ========================================================================
    std::map<FeatureType, int> featureCount = countInstancesByFeature(instances);

    NeighborhoodMgr neighbor_mgr;
    neighbor_mgr.buildFromPairs(neighborPairs, featureCount);

    NRTree orderedNRTree;
    orderedNRTree.build(neighbor_mgr, featureCount,instances);

    // ========================================================================
    // Step 5: Mine Colocation Patterns
    // ========================================================================
    JoinlessMiner miner;

    // Callback đơn giản hơn, không dùng \r để tránh mất log debug
    auto progressCallback = [](int currentStep, int totalSteps, const std::string& message, double percentage) {
        };

    auto colocations = miner.mineColocations(config.minPrev, orderedNRTree, instances, featureCount, progressCallback);

    // ========================================================================
    // Final Report
    // ========================================================================
    auto programEnd = std::chrono::high_resolution_clock::now();
    double totalTime = std::chrono::duration<double>(programEnd - programStart).count();

    // --- REPORT GENERATION (FILE ONLY) ---
    // 1. Get Memory Info (Peak)
    HANDLE handle = GetCurrentProcess();
    PROCESS_MEMORY_COUNTERS memCounter;
    SIZE_T peakMemMB = 0;

    if (GetProcessMemoryInfo(handle, &memCounter, sizeof(memCounter))) {
        peakMemMB = memCounter.PeakWorkingSetSize / 1024 / 1024; // Convert to MB
    }

    // 2. Write to File
    std::ofstream outFile("../results.txt");
    if (!outFile.is_open()) {
        std::cerr << "Cannot open results.txt for writing.\n";
        return 1;
    }
    // (A) Thông tin Dataset & Config
    outFile << "=== FINAL REPORT ===\n";
    outFile << "Dataset Path:      " << config.datasetPath << "\n";
    outFile << "Total Instances:   " << instances.size() << "\n";
    outFile << "Neighbor Distance: " << config.neighborDistance << "\n";
    outFile << "Min Prevalence:    " << config.minPrev << "\n";
    outFile << "----------------------------------------\n";

    // (B) Execution Time
    outFile << "Execution Time: " << std::fixed << std::setprecision(3) << totalTime << " s\n";

    // (C) Peak Memory Usage
    outFile << "Peak Memory Usage: " << peakMemMB << " MB\n";

    // (D) Number of Patterns Found
    outFile << "Patterns Found: " << colocations.size() << "\n";
    outFile << "----------------------------------------\n";

    // (E) List of Patterns
    if (!colocations.empty()) {
        int idx = 1;
        for (const auto& col : colocations) {
            outFile << "[" << idx++ << "] {";
            for (size_t i = 0; i < col.size(); ++i) {
                outFile << (i > 0 ? ", " : "") << col[i];
            }
            outFile << "}\n";
        }
    }
    else {
        outFile << "No patterns found.\n";
    }

    outFile.close();

    std::cout << "Done! Please check 'result.txt'.\n";
    return 0;
}