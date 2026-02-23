/**
 * @file utils.cpp
 * @brief Implementation of utility helper functions
 */

#include "utils.h"
#include "constants.h"
#include <set>
#include <numeric>
#include <chrono>
#include <windows.h>
#include <psapi.h>
#include <iostream> 
#include <iomanip>
#include <algorithm> 
#include <cmath>    

// Get all unique feature types from instances
std::vector<FeatureType> getAllObjectTypes(const std::vector<SpatialInstance>& instances) {
    // Use a set to automatically handle uniqueness
    std::set<FeatureType> objectTypesSet;
    
    // Use std::transform to extract types
    std::transform(instances.begin(), instances.end(),
                   std::inserter(objectTypesSet, objectTypesSet.end()),
                   [](const SpatialInstance& instance) { return instance.type; });
    
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
        FeatureType featureType = instance.id.substr(0, 1);
        featureCount[featureType]++;
    }
    
    return featureCount;
}




std::optional<SpatialInstance> getInstanceByID(
    const std::vector<SpatialInstance>& instances, 
    const instanceID& id) 
{
    // Use std::find_if for search
    const auto it = std::find_if(instances.begin(), instances.end(),
                                  [&id](const SpatialInstance& instance) {
                                      return instance.id == id;
                                  });
    
    // Return instance if found, nullopt otherwise
    return (it != instances.end()) ? std::optional<SpatialInstance>(*it) : std::nullopt;
}

// Step 2: Sorting features in ascending order of the quantity of instances
std::vector<FeatureType> featureSort(const std::vector<FeatureType>& featureSet, const std::vector<SpatialInstance>& instances) {
    // Generate feature counts using the helper function
    const std::map<FeatureType, int> featureCounts = countInstancesByFeature(instances);
    
    // Create a copy to sort (input is const)
    std::vector<FeatureType> sortedFeatures = featureSet;

    // Sort logic based on Algorithm 1 Step 2
    // Ascending order of instance counts
    std::sort(sortedFeatures.begin(), sortedFeatures.end(), 
        [&featureCounts](const FeatureType& a, const FeatureType& b) {
            const int countA = (featureCounts.find(a) != featureCounts.end()) ? featureCounts.at(a) : 0;
            const int countB = (featureCounts.find(b) != featureCounts.end()) ? featureCounts.at(b) : 0;
            
            // Primary sort key: count (ascending)
            if (countA != countB) {
                return countA < countB;
            }
            // Secondary sort key: lexicographical (for stability/determinism)
            return a < b;
        }
    );
    return sortedFeatures;
}

// Step 3: Calculating delta for the spatial dataset
// Formula: delta = (2 / (m*(m-1))) * Sum_{i<j} (num(f_j) / num(f_i))
// This represents the average ratio of instance counts between all pairs of features,
// where features are sorted by instance count (f_i <= f_j).
double calculateDelta(
    const std::vector<FeatureType>& sortedFeatures,
    const std::map<FeatureType, int>& featureCounts) {

    if (sortedFeatures.size() < 2) return 0.0;

    std::vector<double> logCounts;
    logCounts.reserve(sortedFeatures.size());

    // Lấy count theo đúng sortedFeatures
    for (const auto& feat : sortedFeatures) {
        auto it = featureCounts.find(feat);
        if (it != featureCounts.end() && it->second > 0) {
            logCounts.push_back(std::log(static_cast<double>(it->second)));
        }
    }

    size_t m = logCounts.size();
    if (m < 2) return 0.0;

    // Tính mean
    double sumLog = std::accumulate(logCounts.begin(), logCounts.end(), 0.0);
    double meanLog = sumLog / m;

    // Tính variance
    double sumSqDiff = 0.0;
    for (double val : logCounts) {
        sumSqDiff += (val - meanLog) * (val - meanLog);
    }

    double variance = sumSqDiff / (m - 1); // sample stddev

    return std::sqrt(variance);
}


// Calculate Participation Ratio (PR)
// PR(fi, C) = (number of distinct instances of fi in T(C)) / (number of instances of fi)
double calculatePR(
    const FeatureType& featureType,
    const Colocation& pattern,
    const std::map<Colocation, std::vector<ColocationInstance>>& tableInstance,
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

	// Look up tableInstance for the given pattern
    const auto it = tableInstance.find(pattern);
    if (it != tableInstance.end()) {
        const std::vector<ColocationInstance>& instancesList = it->second;

		// Iterate through each row in the table instance
        for (const auto& row : instancesList) {
            if (featureIndex < static_cast<int>(row.size()) && row[featureIndex]) {
                distinctInstances.insert(row[featureIndex]->id);
            }
        }
    }

    // 3. Get total count of featureType globally for denominator
    const int totalCount = (featureCounts.find(featureType) != featureCounts.end()) 
        ? featureCounts.at(featureType) : 0;

    if (totalCount == 0) {
        return 0.0;
    }

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
    const double delta)
{
    if (pattern.empty()) return 0.0;

    // rareType phải thuộc pattern
    if (std::find(pattern.begin(), pattern.end(), rareType) == pattern.end()) {
        return 0.0;
    }

    // 1. Tìm minCount = N(f_min)
    int minCount = -1;
    for (const auto& f : pattern) {
        auto it = featureCounts.find(f);
        if (it != featureCounts.end()) {
            int count = it->second;
            if (minCount == -1 || count < minCount) {
                minCount = count;
            }
        }
    }

    if (minCount <= 0) return 0.0;

    // 2. Lấy count của rareType
    auto itRare = featureCounts.find(rareType);
    if (itRare == featureCounts.end() || itRare->second <= 0) return 0.0;

    int count = itRare->second;

    // 3. Tính theo log-space (GIỐNG calcRareIntensity)
    double sigmaSq2 = 2.0 * delta * delta;
    if (sigmaSq2 <= 0) sigmaSq2 = 1e-9;

    double logMin = std::log(static_cast<double>(minCount));
    double logCount = std::log(static_cast<double>(count));
    double deltaLog = logCount - logMin;

    double ri = std::exp(-(deltaLog * deltaLog) / sigmaSq2);

    return ri;
}


// Calculate Participation Index (PI)
// PI(C) = min_{i=1 to k} { PR(fi, C) }
double calculatePI(
    const Colocation& pattern,
    const std::map<Colocation, std::vector<ColocationInstance>>& tableInstance,
    const std::map<FeatureType, int>& featureCounts) 
{
    if (pattern.empty()) {
        return 0.0;
    }

    double minPR = 1.0; // PR is a probability/ratio <= 1.0
    bool isFirstFeature = true;

    for (const auto& feature : pattern) {
        double pr = calculatePR(feature, pattern, tableInstance, featureCounts);
        if (isFirstFeature) {
            minPR = pr;
            isFirstFeature = false;
        } else {
            if (pr < minPR) {
                minPR = pr;
            }
        }
    }
    
    return minPR;
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
    std::cout << "[PERF] " << stepName << ": " << duration << " ms\n";
}


double getMemoryUsageMB() {
    PROCESS_MEMORY_COUNTERS pmc;
    if (GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc))) {
        return static_cast<double>(pmc.WorkingSetSize) / (1024.0 * 1024.0);
    }
    return 0.0;
}