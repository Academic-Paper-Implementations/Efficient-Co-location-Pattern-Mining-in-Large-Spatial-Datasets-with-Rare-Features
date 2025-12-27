/**
 * @file miner.h
 * @brief Joinless colocation pattern mining algorithm implementation
 * 
 * This file contains the core mining algorithm that discovers prevalent
 * spatial colocation patterns without using expensive join operations.
 */

#pragma once
#include "types.h"
#include "neighborhood_mgr.h"
#include <vector>
#include <map>
#include <functional>

/**
 * @brief Progress callback function type
 * 
 * Callback signature: void(currentStep, totalSteps, message, percentage)
 * Used to report mining progress to the caller.
 */
using ProgressCallback = std::function<void(int, int, const std::string&, double)>;

/**
 * @brief JoinlessMiner class implementing the joinless colocation mining algorithm
 * 
 * This class implements the joinless approach for mining spatial colocation patterns.
 * The algorithm uses star neighborhoods to avoid expensive join operations while
 * discovering prevalent patterns that meet the minimum prevalence threshold.
 */
class JoinlessMiner {
private:
    double minPrev;                          ///< Minimum prevalence threshold
    NRTree* orderedNRTree;        ///< Pointer to neighborhood manager
    ProgressCallback progressCallback;        ///< Progress reporting callback

    std::map<Colocation, std::vector<ColocationInstance>> genTableInstance(
        const std::vector<Colocation>& candidates,
        const std::vector<ColocationInstance>& prevTableInstances,
		const NRTree& orderedNRTree
    );

    // Helper function to find neighbors of an instance for a specific feature type from NRTree
    std::vector<const SpatialInstance*> findNeighbors(
        const NRTree& tree,
        const SpatialInstance* instance,
        const FeatureType& featureType
    );

    // Helper function to calculate S(I, f) = Neigh(o1, f) ∩ ··· ∩ Neigh(ok, f) (Definition 8)
    std::vector<const SpatialInstance*> findExtendedSet(
        const NRTree& tree,
        const ColocationInstance& instance,
        const FeatureType& featureType
    );

    std::vector<ColocationInstance> calculateWPI(
        const std::vector<Colocation>& candidates,
		const std::vector<ColocationInstance>& tableInstances
    );

   
    std::vector<Colocation> selectPrevColocations(
        const std::vector<Colocation>& candidates,
        const std::vector<ColocationInstance>& tableInstances,
        double minPrev,
        const std::map<FeatureType, int>& featureCount, // Cần thêm để tính PR/RI
        double delta
    );




public:
    /**
     * @brief Mine prevalent colocation patterns using the joinless algorithm
     * 
     * Main entry point for the mining algorithm. Discovers all prevalent
     * colocation patterns that meet the minimum prevalence threshold.
     * 
     * @param minPrevalence Minimum prevalence threshold (0.0 to 1.0)
     * @param nbrMgr Pointer to neighborhood manager containing star neighborhoods
     * @param instances Vector of all spatial instances
     * @param progressCb Optional callback for progress reporting
     * @return std::vector<Colocation> All discovered prevalent colocation patterns
     */
    std::vector<Colocation> mineColocations(
        double minPrevalence, 
        NRTree& orderedNRTree, 
        const std::vector<SpatialInstance>& instances,
		const std::map<FeatureType, int>& featureCount,
        ProgressCallback progressCb = nullptr
    );
    
    /**
     * @brief Generate (k+1)-size candidate patterns from k-size prevalent patterns
     * 
     * Uses Apriori-gen approach: joins patterns with matching (k-1) prefixes
     * and prunes candidates whose subsets are not all prevalent.
     * 
     * @param prevPrevalent Vector of k-size prevalent patterns
     * @return std::vector<Colocation> Generated (k+1)-size candidate patterns
     */
    std::vector<Colocation> generateCandidates(
        const std::vector<Colocation>& prevPrevalent,
        const std::map<FeatureType, int>& featureCount
    );

    std::vector<Colocation> filterCandidates(
        const std::vector<Colocation>& candidates,
		const std::vector<Colocation>& prevPrevalent,
		const std::vector<ColocationInstance>& tableInstance,
		double minPrev,
		std::map<FeatureType, int> featureCount,
		double delta
    );


    std::vector<const SpatialInstance*> JoinlessMiner::findNeighbors(
        const NRTree& tree,
        const SpatialInstance* instance,
        const FeatureType& featureType
    );
    std::vector<const SpatialInstance*> JoinlessMiner::findExtendedSet(
        const NRTree& tree,
        const ColocationInstance& instance,
        const FeatureType& featureType
    );
};