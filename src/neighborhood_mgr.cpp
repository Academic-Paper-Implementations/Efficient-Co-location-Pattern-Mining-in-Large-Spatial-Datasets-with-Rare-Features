/**
 * @file neighborhood_mgr.cpp
 * @brief Implementation of ordered neighborhood management
 */

#include "neighborhood_mgr.h"
#include "utils.h"
#include <algorithm>


/**
 * @brief Check if neighbor should be included in center's ordered neighborhood
 * @param centerType Feature type of the center instance
 * @param neighborType Feature type of the neighbor instance
 * @param counts Map of instance counts for each feature type
 * @return bool True if neighbor should be in center's ordered neighborhood
 * 
 * Ordering is based on instance count (ascending). For equal counts, uses lexicographic order.
 */
bool NeighborhoodMgr::isOrdered(const FeatureType& centerType,
    const FeatureType& neighborType,
    const std::map<FeatureType, int>& counts){

    int numCenter = counts.at(centerType);
    int numNeighbor = counts.at(neighborType);

    if (numCenter < numNeighbor) return true;
    if (numCenter == numNeighbor) return centerType <= neighborType;
    return false;
}


/**
 * @brief Build ordered neighborhoods from neighbor pairs
 * @param pairs Vector of neighbor pairs found by spatial indexing
 * @param featureCounts Map of instance counts for each feature type
 * 
 * Constructs bidirectional ordered neighborhoods. For each pair (A, B):
 * - If A's feature count <= B's feature count, B is added to A's neighborhood
 * - If B's feature count <= A's feature count, A is added to B's neighborhood
 */
void NeighborhoodMgr::buildFromPairs(const std::vector<std::pair<SpatialInstance, SpatialInstance>>& pairs,
    const std::map<FeatureType, int>& featureCounts) {
    
    orderedNeighborMap.clear();
    
    for (const auto& pair : pairs) {
        const SpatialInstance& center = pair.first;
        const SpatialInstance& neighbor = pair.second;

        // Check if neighbor belongs to center's ordered neighborhood
        if (isOrdered(center.type, neighbor.type, featureCounts)) {
            auto& neighborhoodList = orderedNeighborMap[center.type];
            auto existingNeighborhood = std::find_if(neighborhoodList.begin(), neighborhoodList.end(), [&](const OrderedNeigh& set) {
                return set.center->id == center.id;
                });
                
            if (existingNeighborhood != neighborhoodList.end()) {
                existingNeighborhood->neighbors[neighbor.type].push_back(&neighbor);
            }
            else {
                OrderedNeigh newSet;
                newSet.center = &center;
                newSet.neighbors[neighbor.type].push_back(&neighbor);
                neighborhoodList.push_back(newSet);
            }
        }
        
        // Check if center belongs to neighbor's ordered neighborhood
        if (isOrdered(neighbor.type, center.type, featureCounts)) {
            auto& neighborhoodList = orderedNeighborMap[neighbor.type];
            auto existingNeighborhood = std::find_if(neighborhoodList.begin(), neighborhoodList.end(), [&](const OrderedNeigh& set) {
                return set.center->id == neighbor.id;
                });
                
            if (existingNeighborhood != neighborhoodList.end()) {
                existingNeighborhood->neighbors[center.type].push_back(&center);
            }
            else {
                OrderedNeigh newSet;
                newSet.center = &neighbor;
                newSet.neighbors[center.type].push_back(&center);
                neighborhoodList.push_back(newSet);
            }
        }
    }
}


/**
 * @brief Get the ordered neighborhood map
 * @return const reference to the ordered neighborhood map
 */
const std::unordered_map<FeatureType, std::vector<OrderedNeigh>>& NeighborhoodMgr::getOrderedNeighbors() const {
    return orderedNeighborMap;
}