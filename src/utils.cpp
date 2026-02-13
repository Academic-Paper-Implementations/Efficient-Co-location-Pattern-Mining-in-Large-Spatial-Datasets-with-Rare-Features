/**
 * @file utils.cpp
 * @brief Implementation of utility helper functions
 */

#include "utils.h"
#include "constants.h"
#include <set>
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
double calculateDelta(const std::vector<FeatureType>& sortedFeatures, const std::map<FeatureType, int>& featureCounts) {
    if (sortedFeatures.size() < 2) {
        return 0.0;
    }

    // 1. Extract counts in the order of sortedFeatures (Step 2 order)
    std::vector<double> counts;
    counts.reserve(sortedFeatures.size());
    for (const auto& feat : sortedFeatures) {
        if (featureCounts.find(feat) != featureCounts.end()) {
            counts.push_back(static_cast<double>(featureCounts.at(feat)));
        } else {
            // Should not happen if sortedFeatures comes from keys of featureCounts,
            // but safe to handle.
            counts.push_back(0.0);
        }
    }

    // 2. NO sorting of 'counts' here. 
    // We rely on 'sortedFeatures' being already sorted by quantity (Step 2).
    // Paper: delta = 2/(m(m-1)) * Sum_{i<j} (|fj| / |fi|)
    // The indices i, j correspond to the sorted feature list order.

    const double numFeatures = static_cast<double>(counts.size());
    double sumRatios = 0.0;

    // 3. Calculate sum of ratios for all pairs i < j
    for (size_t i = 0; i < counts.size(); ++i) {
        for (size_t j = i + 1; j < counts.size(); ++j) {
            // Formula uses |fj| / |fi| where i < j
            // Since features are sorted by count ascending, |fi| <= |fj| usually holds,
            // making the ratio >= 1 (or close to it/handling stability).
            const double numerator = counts[j];
            double denominator = counts[i];
            
            // Handle division by zero
            if (denominator == 0.0) {
                denominator = Constants::EPSILON_SMALL;
            }
            
            const double ratio = numerator / denominator;
            sumRatios += ratio;
        }
    }

    // 4. Calculate final Delta
    // Factor = 2 / (numFeatures * (numFeatures - 1))
    const double factor = 2.0 / (numFeatures * (numFeatures - 1.0));
    
    return factor * sumRatios;
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
    // Safety check for delta to avoid division by zero
    if (delta <= Constants::EPSILON_DELTA) return 0.0;

    // Definition check: RI(fi, C) is only defined if fi is in C
    if (std::find(pattern.begin(), pattern.end(), rareType) == pattern.end()) {
        return 0.0;
    }

    // 1. Find num(f_min) in the pattern
    int minCount = -1;
    for (const auto& feature : pattern) {
        const auto featureIt = featureCounts.find(feature);
        if (featureIt != featureCounts.end()) {
            const int count = featureIt->second;
            if (minCount == -1 || count < minCount) {
                minCount = count;
            }
        } else {
            // If a feature in the pattern has 0 instances, minCount is 0 logic
            minCount = 0; 
            break; 
        }
    }

    if (minCount <= 0) {
        return 0.0; // Avoid division by zero in v calculation
    }

    // 2. Get num(fi) for the rareType
    const int rareCount = (featureCounts.find(rareType) != featureCounts.end()) 
        ? featureCounts.at(rareType) : 0;

    // 3. Calculate v(fi, C) = num(fi) / num(f_min)
    const double v_val = static_cast<double>(rareCount) / static_cast<double>(minCount);

    // 4. Calculate RI
    // exponent term = - (v - 1)^2 / (2 * delta^2)
    const double numerator = std::pow(v_val - 1.0, 2);
    const double denominator = 2.0 * delta * delta;
    
    return std::exp(-numerator / denominator);
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