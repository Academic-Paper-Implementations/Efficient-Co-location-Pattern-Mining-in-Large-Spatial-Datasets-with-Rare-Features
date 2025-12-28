/**
 * @file spatial_index.cpp
 * @brief Implementation of spatial indexing and neighbor search
 */

#include "spatial_index.h"
#include <cmath>
#include <iostream>
#include <algorithm>


/**
 * @brief Constructor to initialize SpatialIndex with a distance threshold
 * @param distThresh Maximum distance for two instances to be considered neighbors
 */
SpatialIndex::SpatialIndex(double distThresh)
    : distanceThreshold(distThresh)
{
}


/**
 * @brief Calculate Euclidean distance between two spatial instances
 * @param a First spatial instance
 * @param b Second spatial instance
 * @return double Euclidean distance between a and b
 */
double SpatialIndex::euclideanDist(const SpatialInstance& a, const SpatialInstance& b) {
    double dx = a.x - b.x;
    double dy = a.y - b.y;
    return std::sqrt(dx * dx + dy * dy);
}


/**
 * @brief Find all neighbor pairs within the distance threshold
 * @param instances Vector of all spatial instances to search
 * @return std::vector<std::pair<SpatialInstance, SpatialInstance>> Vector of neighbor pairs
 * 
 * Uses grid-based spatial partitioning to optimize neighbor search from O(n²) to O(n).
 * Divides the spatial domain into grid cells and only checks instances in adjacent cells.
 */
std::vector<std::pair<SpatialInstance, SpatialInstance>> SpatialIndex::findNeighborPair(const std::vector<SpatialInstance>& instances) {
    std::vector<std::pair<SpatialInstance, SpatialInstance>> neighborPairs;

    // Calculate spatial bounds
    double minX = std::min_element(instances.begin(), instances.end(),
        [](const SpatialInstance& a, const SpatialInstance& b) { return a.x < b.x; })->x;
    double minY = std::min_element(instances.begin(), instances.end(),
        [](const SpatialInstance& a, const SpatialInstance& b) { return a.y < b.y; })->y;
    double maxX = std::max_element(instances.begin(), instances.end(),
        [](const SpatialInstance& a, const SpatialInstance& b) { return a.x < b.x; })->x;
    double maxY = std::max_element(instances.begin(), instances.end(),
        [](const SpatialInstance& a, const SpatialInstance& b) { return a.y < b.y; })->y;

    // Create grid cells based on distance threshold
    int gridX = static_cast<int>(std::ceil((maxX - minX) / distanceThreshold));
    int gridY = static_cast<int>(std::ceil((maxY - minY) / distanceThreshold));
    std::vector<std::vector<SpatialInstance>> gridCells(gridX * gridY);

    // Assign instances to grid cells
    for (const auto& inst : instances) {
        int cx = static_cast<int>((inst.x - minX) / distanceThreshold);
        int cy = static_cast<int>((inst.y - minY) / distanceThreshold);
        gridCells[cx * gridY + cy].push_back(inst);
    }

    // Check pairs within and between adjacent cells
    for (int cx = 0; cx < gridX; ++cx) {
        for (int cy = 0; cy < gridY; ++cy) {
            const auto& cell = gridCells[cx * gridY + cy];

            // Check pairs within the same cell
            for (size_t i = 0; i < cell.size(); ++i) {
                for (size_t j = i + 1; j < cell.size(); ++j) {
                    if (cell[i].type != cell[j].type && euclideanDist(cell[i], cell[j]) <= distanceThreshold) {
                        neighborPairs.emplace_back(cell[i], cell[j]);
                    }
                }
                
                // Check pairs with adjacent cells (avoid duplicate checks)
                for (int dx = 0; dx <= 1; ++dx) {
                    for (int dy = (dx == 0 ? 1 : -1); dy <= 1; ++dy) {
                        if (dx == 0 && dy == 0) continue;
                        
                        int nx = cx + dx;
                        int ny = cy + dy;
                        
                        if (nx >= 0 && nx < gridX && ny >= 0 && ny < gridY) {
                            const auto& neighborCell = gridCells[nx * gridY + ny];
                            for (const auto& neighborInst : neighborCell) {
                                if (cell[i].type != neighborInst.type && euclideanDist(cell[i], neighborInst) <= distanceThreshold) {
                                    neighborPairs.emplace_back(cell[i], neighborInst);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    return neighborPairs;
}