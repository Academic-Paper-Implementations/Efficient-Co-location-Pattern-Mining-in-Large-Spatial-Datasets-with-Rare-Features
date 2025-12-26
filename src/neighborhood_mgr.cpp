/**
 * @file neighborhood_mgr.cpp
 * @brief Implementation of star neighborhood management
 */

#include "neighborhood_mgr.h"
#include "utils.h"
#include <algorithm>

//Check is ordered function
bool NeighborhoodMgr::isOrdered(const FeatureType& centerType,
    const FeatureType& neighborType,
    const std::map<FeatureType, int>& counts){

    int numCenter = counts.at(centerType);
    int numNeighbor = counts.at(neighborType);

    if (numCenter < numNeighbor) return true;
    if (numCenter == numNeighbor) return centerType <= neighborType;
    return false;
}

void NeighborhoodMgr::buildFromPairs(const std::vector<std::pair<SpatialInstance, SpatialInstance>>& pairs,
    const std::map<FeatureType, int>& featureCounts) {
    // Build star neighborhoods from neighbor pairs
    // A star neighborhood has a center instance and all its neighbors
    orderedNeighborMap.clear();
    for (const auto& pair : pairs) {
        const SpatialInstance& center = pair.first;
        const SpatialInstance& neighbor = pair.second;

        // --- CHIỀU 1: Kiểm tra xem neighbor có thuộc OrderedNeighborset của center không? ---
        if (isOrdered(center.type, neighbor.type, featureCounts)) {
            auto& vec = orderedNeighborMap[center.type];
            auto it = std::find_if(vec.begin(), vec.end(), [&](const OrderedNeigh& set) {
                return set.center->id == center.id;
                });
            if (it != vec.end()) {
                it->neighbors[neighbor.type].push_back(&neighbor);
            }
            else {
                OrderedNeigh newSet;
                newSet.center = &center;
                newSet.neighbors[neighbor.type].push_back(&neighbor);
                vec.push_back(newSet);
            }
        }
        // --- CHIỀU 2: Kiểm tra xem center có thuộc OrderedNeighborset của neighbor không? ---
        if (isOrdered(neighbor.type, center.type, featureCounts)) {
            auto& vec = orderedNeighborMap[neighbor.type];
            auto it = std::find_if(vec.begin(), vec.end(), [&](const OrderedNeigh& set) {
                return set.center->id == neighbor.id;
                });
            if (it != vec.end()) {
                it->neighbors[center.type].push_back(&center);
            }
            else {
                OrderedNeigh newSet;
                newSet.center = &neighbor;
                newSet.neighbors[center.type].push_back(&center);
                vec.push_back(newSet);
            }
        }
    }
}
const std::unordered_map<FeatureType, std::vector<OrderedNeigh>>& NeighborhoodMgr::getOrderedNeighbors() const {
    return orderedNeighborMap;
}