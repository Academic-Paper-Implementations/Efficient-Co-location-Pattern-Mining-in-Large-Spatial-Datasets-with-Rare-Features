/**
 * @file utils.cpp
 * @brief Implementation of utility helper functions
 */

#include "utils.h"
#include <set>
#include <chrono>
#include <windows.h>
#include <psapi.h>
#include <iostream> 
#include <iomanip>
#include <windows.h>
#include <psapi.h>
#include <algorithm> 
#include <cmath>    

// Get all unique feature types from instances
std::vector<FeatureType> getAllObjectTypes(const std::vector<SpatialInstance>& instances) {
    // Use a set to automatically handle uniqueness
    std::set<FeatureType> objectTypesSet;
    
    for (const auto& instance : instances) {
        FeatureType objectType = instance.type;
        objectTypesSet.insert(objectType);
    }
    
    // Convert set to vector (set maintains sorted order)
    return std::vector<FeatureType>(objectTypesSet.begin(), objectTypesSet.end());
}

// Count the number of instances for each feature type
std::map<FeatureType, int> countInstancesByFeature(const std::vector<SpatialInstance>& instances) {
    std::map<FeatureType, int> featureCount;
    
    for (const auto& instance : instances) {
        // Extract feature type from instance ID
        // Assumes ID format is: FeatureType + Number (e.g., "A1", "B2")
        // Takes first character as feature type
        FeatureType objectType = instance.id.substr(0, 1);
        featureCount[objectType]++;
    }
    
    return featureCount;
}


SpatialInstance getInstanceByID(
    const std::vector<SpatialInstance>& instances, 
    const instanceID& id) 
{
    // Linear search for instance with matching ID
    for (const auto& instance : instances) {
        if (instance.id == id) {
            return instance;
        }
    }
    
    // Return empty instance if not found
    return SpatialInstance{};
}

// Step 2: Sorting features in ascending order of the quantity of instances
std::vector<FeatureType> featureSort(std::vector<FeatureType>& featureSet, const std::vector<SpatialInstance>& instances) {
    // Generate feature counts using the helper function
    std::map<FeatureType, int> featureCounts = countInstancesByFeature(instances);

    // Sort logic based on Algorithm 1 Step 2
    // Ascending order of instance counts
    std::sort(featureSet.begin(), featureSet.end(), 
        [&featureCounts](const FeatureType& a, const FeatureType& b) {
            int countA = 0;
            int countB = 0;
            if (featureCounts.find(a) != featureCounts.end()) countA = featureCounts.at(a);
            if (featureCounts.find(b) != featureCounts.end()) countB = featureCounts.at(b);
            
            // Primary sort key: count (ascending)
            if (countA != countB) {
                return countA < countB;
            }
            // Secondary sort key: lexicographical (for stability/determinism)
            return a < b;
        }
    );
    return featureSet;
}

// Step 3: Calculating delta for the spatial dataset
// Formula: delta = (2 / (m*(m-1))) * Sum_{i<j} (num(f_j) / num(f_i))
// This represents the average ratio of instance counts between all pairs of features,
// where features are sorted by instance count (f_i <= f_j).
double calculateDelta(const std::map<FeatureType, int>& featureCounts) {
    if (featureCounts.size() < 2) {
        return 0.0;
    }

    // 1. Extract counts into a vector for indexed access
    std::vector<double> counts;
    for (const auto& pair : featureCounts) {
        counts.push_back(static_cast<double>(pair.second));
    }

    // 2. Sort counts in ascending order (matching Step 2's logic implication)
    // This ensures that for i < j, counts[i] <= counts[j], so the ratio >= 1.
    std::sort(counts.begin(), counts.end());

    double m = static_cast<double>(counts.size());
    double sumRatios = 0.0;

    // 3. Calculate sum of ratios for all pairs i < j
    for (int i = 0; i < m; ++i) {
        for (int j = i + 1; j < m; ++j) {
            // Avoid division by zero if count is 0 (though rare in mining)
            double denominator = counts[i];
            if (denominator == 0) denominator = 1e-6; 
            
            double ratio = counts[j] / denominator;
            sumRatios += ratio;
        }
    }

    // 4. Calculate final Delta
    // Factor = 2 / (m * (m - 1))
    double factor = 2.0 / (m * (m - 1.0));
    
    return factor * sumRatios;
}

// Calculate Participation Ratio (PR)
// PR(fi, C) = (number of distinct instances of fi in T(C)) / (number of instances of fi)
double calculatePR(
    const FeatureType& featureType,
    const Colocation& pattern,
    const std::vector<ColocationInstance>& tableInstance,
    const std::map<FeatureType, int>& featureCounts) 
{
    // 1. Find the index of featureType in the pattern
    int featureIndex = -1;
    for (size_t i = 0; i < pattern.size(); ++i) {
        if (pattern[i] == featureType) {
            featureIndex = static_cast<int>(i);
            break;
        }
    }

    if (featureIndex == -1) {
        // Feature not in pattern
        return 0.0;
    }

    // 2. Count distinct instances of featureType in T(C) to get numerator
    std::set<instanceID> distinctInstances;
    for (const auto& row : tableInstance) {
        // Safety check: row size should match pattern size
        if (featureIndex < static_cast<int>(row.size()) && row[featureIndex]) {
            distinctInstances.insert(row[featureIndex]->id);
        }
    }

    // 3. Get total count of featureType globally for denominator
    int totalCount = 0;
    if (featureCounts.find(featureType) != featureCounts.end()) {
        totalCount = featureCounts.at(featureType);
    }

    if (totalCount == 0) return 0.0;

    // 4. Calculate PR
    return static_cast<double>(distinctInstances.size()) / static_cast<double>(totalCount);
}

// Calculate Rare Intensity (RI) for a feature in a co-location pattern
// Definition 3, Formula (5): RI(fi, C) = exp( - (v(fi, C) - 1)^2 / (2 * delta^2) )
// where v(fi, C) = num(fi) / num(f_min) (Definition 2)
double calculateRareIntensity(
    const FeatureType& rareType, 
    const Colocation& pattern,
    const std::map<FeatureType, int>& featureCounts,
    double delta) 
{
    // Safety check for delta to avoid division by zero
    if (delta <= 1e-9) return 0.0;

    // 1. Find num(f_min) in the pattern
    int minCount = -1;
    for (const auto& feature : pattern) {
        if (featureCounts.find(feature) != featureCounts.end()) {
            int count = featureCounts.at(feature);
            if (minCount == -1 || count < minCount) {
                minCount = count;
            }
        } else {
            // If a feature in the pattern has 0 instances, minCount is 0 logic
            minCount = 0; 
            break; 
        }
    }

    if (minCount <= 0) return 0.0; // Avoid division by zero in v calculation

    // 2. Get num(fi) for the rareType
    int rareCount = 0;
    if (featureCounts.find(rareType) != featureCounts.end()) {
        rareCount = featureCounts.at(rareType);
    }

    // 3. Calculate v(fi, C) = num(fi) / num(f_min)
    double v_val = static_cast<double>(rareCount) / static_cast<double>(minCount);

    // 4. Calculate RI
    // exponent term = - (v - 1)^2 / (2 * delta^2)
    double numerator = std::pow(v_val - 1.0, 2);
    double denominator = 2.0 * delta * delta;
    
    return std::exp(-numerator / denominator);
}
void findCombinations(
    const std::vector<FeatureType>& candidatePattern,
    int typeIndex,
    std::vector<const SpatialInstance*>& currentInstance,
    const std::unordered_map<FeatureType, std::vector<const SpatialInstance*>>& neighborMap,
    std::vector<ColocationInstance>& results) 
{
    // Base case: if we've matched all types in the candidate pattern
    if (typeIndex >= candidatePattern.size()) {
        results.push_back(currentInstance);
        return;
    }
    FeatureType currentType = candidatePattern[typeIndex];

    auto it = neighborMap.find(currentType);
    if (it != neighborMap.end()) {
        for (const auto* neighbor : it->second) {
            currentInstance.push_back(neighbor);
            findCombinations(candidatePattern, typeIndex + 1, currentInstance, neighborMap, results);
            currentInstance.pop_back();
        }
    }
}


void printDuration(const std::string& stepName, std::chrono::high_resolution_clock::time_point start, std::chrono::high_resolution_clock::time_point end) {
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "[PERF] " << stepName << ": " << duration << " ms" << std::endl;
}


double getMemoryUsageMB() {
    PROCESS_MEMORY_COUNTERS pmc;
    if (GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc))) {
        return static_cast<double>(pmc.WorkingSetSize) / (1024.0 * 1024.0);
    }
    return 0.0;
}