/**
 * @file utils.h
 * @brief Utility helper functions for spatial colocation mining
 */

#pragma once
#include "types.h"
#include <vector>
#include <set>
#include <string>
#include <chrono>
#include <map>

/**
 * @brief Get all unique feature types from spatial instances
 * 
 * Extracts and returns a sorted vector of all unique feature types
 * present in the given instances.
 * 
 * @param instances Vector of spatial instances
 * @return std::vector<FeatureType> Sorted vector of unique feature types
 */
std::vector<FeatureType> getAllObjectTypes(const std::vector<SpatialInstance>& instances);

/**
 * @brief Count the number of instances for each feature type
 * 
 * Creates a frequency map showing how many instances exist for each feature type.
 * 
 * @param instances Vector of spatial instances
 * @return std::map<FeatureType, int> Map from feature type to instance count
 */
std::map<FeatureType, int> countInstancesByFeature(const std::vector<SpatialInstance>& instances);

/**
 * @brief Find a spatial instance by its ID
 * 
 * Searches for an instance with the specified ID. Returns an empty
 * SpatialInstance struct if no matching instance is found.
 * 
 * @param instances Vector of spatial instances to search
 * @param id Instance ID to find
 * @return SpatialInstance The found instance, or empty struct if not found
 */
SpatialInstance getInstanceByID(
    const std::vector<SpatialInstance>& instances, 
    const instanceID& id);

/**
 * @brief Sort features by instance count (ascending)
 * @param featureSet The feature set to sort
 * @param instances The instances used to count feature frequency
 */
std::vector<FeatureType> featureSort(std::vector<FeatureType>& featureSet, const std::vector<SpatialInstance>& instances);

double calculateDelta(const std::vector<FeatureType>& sortedFeatures, const std::map<FeatureType, int>& featureCounts);

/**
 * @brief Calculate Participation Ratio (PR) for a feature in a co-location
 * 
 * @param featureType The feature type to calculate PR for
 * @param pattern The co-location pattern
 * @param tableInstance The table instance T(C) of the pattern
 * @param featureCounts Map of total instance counts for all features
 * @return double Participation Ratio value
 */
double calculatePR(
    const FeatureType& featureType,
    const Colocation& pattern,
    const std::vector<ColocationInstance>& tableInstance,
    const std::map<FeatureType, int>& featureCounts);
/**
 * @brief Calculate Rare Intensity (RI) for a feature in a pattern
 * 
 * @param rareType The feature type to calculate RI for
 * @param pattern The co-location pattern
 * @param featureCounts Map of instance counts for all features
 * @param delta Global degree of dispersion
 * @return double Rare Intensity value
 */
double calculateRareIntensity(
    const FeatureType& rareType, 
    const Colocation& pattern,
    const std::map<FeatureType, int>& featureCounts,
    double delta);
/**
* @brief Recursive helper to find all combinations of spatial instances
*        matching a candidate pattern within a star neighborhood.
* 
* @param candidatePattern The candidate colocation pattern being matched
* @param typeIndex Current index in the candidate pattern being processed
* @param currentInstance Current partial instance being built
* @param neighborMap Map of feature types to their neighboring instances
* @param results Vector to store the resulting colocation instances
*/
void findCombinations(
    const std::vector<FeatureType>& candidatePattern,
    int typeIndex,
    std::vector<const SpatialInstance*>& currentInstance,
    const std::unordered_map<FeatureType, std::vector<const SpatialInstance*>>& neighborMap,
    std::vector<ColocationInstance>& results);


/**
* @brief Print the duration of a processing step
* 
* @param stepName Name of the processing step
* @param start Start time point
* @param end End time point
*/
void printDuration(const std::string& stepName, std::chrono::high_resolution_clock::time_point start, std::chrono::high_resolution_clock::time_point end);


/**
 * @brief Get current memory usage in megabytes
 *
 * Uses Windows API to retrieve the current process's memory usage.
 *
 * @return double Memory usage in MB
 */
double getMemoryUsageMB();