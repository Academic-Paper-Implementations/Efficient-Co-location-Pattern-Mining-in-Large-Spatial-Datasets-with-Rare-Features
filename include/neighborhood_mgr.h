/**
 * @file neighborhood_mgr.h
 * @brief Star neighborhood management for spatial instances
 */

#pragma once
#include "types.h"
#include <unordered_map>
#include <vector>
#include "NRTree.h"

/**
 * @brief NeighborhoodMgr class for managing star neighborhoods of spatial instances
 * 
 * Organizes spatial instances into star neighborhoods, where each star consists of
 * a center instance and all its neighbors within the distance threshold.
 */
class NeighborhoodMgr {
public:
    /**
     * @brief Build star neighborhoods from neighbor pairs
     * 
     * Constructs star neighborhoods by grouping neighbor pairs. For each instance,
     * creates a star with that instance as center and all its neighbors.
     * 
     * @param pairs Vector of neighbor pairs found by spatial indexing
     */
    void buildFromPairs(const std::vector<std::pair<SpatialInstance, SpatialInstance>>& pairs,
                        const std::unordered_map<FeatureType, int>&FeatureCounts);
    
    /**
     * @brief Get all ordered neighborhoods organized by feature type
     * 
     * @return const std::unordered_map<FeatureType, std::vector<OrderedNeigh>>& 
     *         Map from feature type to vector of ordered neighborhoods
     */
    const std::unordered_map<FeatureType, std::vector<OrderedNeigh>>& getOrderedNeighbors() const;

private:
    /**
     * @brief Lưu trữ Neigh: Map từ FeatureType -> Danh sách các OrderedNeigh
     */
    std::unordered_map<FeatureType, std::vector<OrderedNeigh>> orderedNeighborMap;

    // Hàm kiểm tra thứ tự
    bool isOrdered(const FeatureType& centerType,
        const FeatureType& neighborType,
        const std::unordered_map<FeatureType, int>& counts);
};