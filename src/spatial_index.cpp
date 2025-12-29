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
double SpatialIndex::euclideanDist(const SpatialInstance& a, const SpatialInstance& b) const {
    const double dx = a.x - b.x;
    const double dy = a.y - b.y;
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
std::vector<std::pair<SpatialInstance, SpatialInstance>> SpatialIndex::findNeighborPair(const std::vector<SpatialInstance>& instances) const {
    std::vector<std::pair<SpatialInstance, SpatialInstance>> neighborPairs;

    // Safety check: empty instances
    if (instances.empty()) {
        return neighborPairs;
    }

    // Calculate spatial bounds
    const double minX = std::min_element(instances.begin(), instances.end(),
        [](const SpatialInstance& a, const SpatialInstance& b) { return a.x < b.x; })->x;
    const double minY = std::min_element(instances.begin(), instances.end(),
        [](const SpatialInstance& a, const SpatialInstance& b) { return a.y < b.y; })->y;
    const double maxX = std::max_element(instances.begin(), instances.end(),
        [](const SpatialInstance& a, const SpatialInstance& b) { return a.x < b.x; })->x;
    const double maxY = std::max_element(instances.begin(), instances.end(),
        [](const SpatialInstance& a, const SpatialInstance& b) { return a.y < b.y; })->y;

    // Create grid cells based on distance threshold
    const int gridCellsX = static_cast<int>(std::ceil((maxX - minX) / distanceThreshold));
    const int gridCellsY = static_cast<int>(std::ceil((maxY - minY) / distanceThreshold));
    std::vector<std::vector<SpatialInstance>> gridCells(gridCellsX * gridCellsY);

    // Assign instances to grid cells
    for (const auto& inst : instances) {
        const int cellX = static_cast<int>((inst.x - minX) / distanceThreshold);
        const int cellY = static_cast<int>((inst.y - minY) / distanceThreshold);
        gridCells[cellX * gridCellsY + cellY].push_back(inst);
    }

    // Check pairs within and between adjacent cells
    for (int cellX = 0; cellX < gridCellsX; ++cellX) {
        for (int cellY = 0; cellY < gridCellsY; ++cellY) {
            const auto& cell = gridCells[cellX * gridCellsY + cellY];

            // Check pairs within the same cell
            for (size_t i = 0; i < cell.size(); ++i) {
                for (size_t j = i + 1; j < cell.size(); ++j) {
                    if (cell[i].type != cell[j].type && euclideanDist(cell[i], cell[j]) <= distanceThreshold) {
                        neighborPairs.emplace_back(cell[i], cell[j]);
                    }
                }
                
                // Check pairs with adjacent cells (avoid duplicate checks)
                for (int deltaX = 0; deltaX <= 1; ++deltaX) {
                    for (int deltaY = (deltaX == 0 ? 1 : -1); deltaY <= 1; ++deltaY) {
                        if (deltaX == 0 && deltaY == 0) {
                            continue;
                        }
                        
                        const int neighborCellX = cellX + deltaX;
                        const int neighborCellY = cellY + deltaY;
                        
                        if (neighborCellX >= 0 && neighborCellX < gridCellsX && neighborCellY >= 0 && neighborCellY < gridCellsY) {
                            const auto& neighborCell = gridCells[neighborCellX * gridCellsY + neighborCellY];
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