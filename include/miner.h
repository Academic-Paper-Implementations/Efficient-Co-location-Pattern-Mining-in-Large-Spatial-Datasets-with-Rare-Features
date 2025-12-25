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

    std::vector<ColocationInstance> genTableInstance(
        const std::vector<Colocation>& candidates,
        const std::vector<ColocationInstance> prevTableInstances,
		const NRTree& orderedNRTree
    );

    std::vector<ColocationInstance> calculateWPI(
        const std::vector<Colocation>& candidates,
		const std::vector<ColocationInstance>& tableInstances
    );

   
    std::vector<Colocation> selectPrevColocations(
        const std::vector<Colocation>& candidates,
        const std::vector<ColocationInstance>& tableInstances,
        double minPrev
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
        const std::vector<Colocation>& prevPrevalent
    );

    std::vector<Colocation> filterCandidates(
        const std::vector<Colocation>& candidates,
		const std::vector<Colocation>& prevPrevalent
	);
};