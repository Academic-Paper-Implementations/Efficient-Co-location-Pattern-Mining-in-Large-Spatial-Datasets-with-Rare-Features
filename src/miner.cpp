/**
 * @file miner.cpp
 * @brief Implementation of the ordered NR-tree based mining algorithm
 * * This file implements the core mining algorithm that discovers prevalent spatial
 * co-location patterns using the Ordered NR-Tree structure as described in Algorithm 1.
 */

#include "miner.h"
#include "utils.h"
#include "neighborhood_mgr.h"
#include "types.h"
#include <algorithm>
#include <set>
#include <string>
#include <iostream>
#include <omp.h> 
#include <iomanip>
#include <chrono>

std::vector<Colocation> JoinlessMiner::mineColocations(
    double minPrev,
    NRTree& orderedNRTree,
    const std::vector<SpatialInstance>& instances,
    const std::map<FeatureType, int>& featureCount,
    ProgressCallback progressCb
) {
    // Start timer
    auto minerStart = std::chrono::high_resolution_clock::now();

    // Assign parameters to member variables for use in other methods
    this->progressCallback = progressCb;

    // Variables: (a) delta; (b) k; (c) Ck; (d) Tk; (e) Pk; (f) Neigh
    // Steps 1-5: Initialization (counting instances, sorting features, calculating delta, gen_Neigh, gen_ordered-NR-tree)
    // Note: Passed in via arguments or pre-calculated in caller
    int k = 2;  // Start with size-2 patterns
    std::vector<FeatureType> types = getAllObjectTypes(instances);
	std::vector<FeatureType> sortedTypes = featureSort(types, instances);
	double delta = calculateDelta(sortedTypes, featureCount);


    // P1 = F (Set of size-1 prevalent co-locations)
    std::vector<Colocation> prevColocations;
    // T1 = O (Table instance of size-1)
    std::vector<ColocationInstance> prevTableInstances;
    std::vector<Colocation> allPrevalentColocations;

    // Estimate total iterations (max pattern size is number of types)
    int maxK = types.size();
    int currentIteration = 0;
    int totalIterations = 0;

    if (progressCallback) {
        progressCallback(0, maxK, "Initializing mining process (Steps 1-6)...", 0.0);
    }

    // Step 6: let P1 = F, T1 = O, k = 2
    for (auto t : types) prevColocations.push_back({ t });

    // Step 7: while Pk-1 not empty do
    while (!prevColocations.empty()) {
        currentIteration++;
        totalIterations = currentIteration;

        // Calculate progress
        double progressPercent = std::min(95.0, (static_cast<double>(currentIteration) / maxK) * 95.0);

        if (progressCallback) {
            progressCallback(currentIteration, maxK,
                "Processing iteration k=" + std::to_string(k) + "...",
                progressPercent);
        }
        std::vector<ColocationInstance> tableInstances;

        // Step 8: Ck = gen_candidate_patterns(Pk-1, k)
        auto t1_start = std::chrono::high_resolution_clock::now();
        std::vector<Colocation> candidates = generateCandidates(prevColocations, featureCount);
        auto t1_end = std::chrono::high_resolution_clock::now();
        printDuration("Step 8: gen_candidate_patterns (k=" + std::to_string(k) + ")", t1_start, t1_end);

        if (candidates.empty()) {
            if (progressCallback) {
                progressCallback(currentIteration, maxK,
                    "No more candidates found. Mining completed.",
                    100.0);
            }
            break;
        }

        if (progressCallback) {
            double progressPercent = std::min(95.0, (static_cast<double>(currentIteration) / maxK) * 95.0);
            progressCallback(currentIteration, maxK,
                "Filtering candidates (Lemma 2 & 3)...",
                progressPercent);
        }

        // Step 9: filter_candidate_patterns(Ck, Pk-1)
        // Uses Lemma 2 and Lemma 3 to prune search space
        auto t2_start = std::chrono::high_resolution_clock::now();
		std::vector<Colocation> fiteredCandidates = candidates;
        if (k==2){
            fiteredCandidates = filterCandidates(candidates, prevColocations, prevTableInstances, minPrev, featureCount, delta);
        }
        auto t2_end = std::chrono::high_resolution_clock::now();
        printDuration("Step 9: filter_candidate_patterns (k=" + std::to_string(k) + ")", t2_start, t2_end);

        // Step 10: Tk = gen_table_instances(Ck, Tk-1, ordered-NR-tree)
        auto t3_start = std::chrono::high_resolution_clock::now();
        tableInstances = genTableInstance(
            fiteredCandidates,
            prevTableInstances,
            orderedNRTree
        );
        auto t3_end = std::chrono::high_resolution_clock::now();
        printDuration("Step 10: gen_table_instances (k=" + std::to_string(k) + ")", t3_start, t3_end);

        if (progressCallback) {
            double progressPercent = std::min(95.0, (static_cast<double>(currentIteration) / maxK) * 95.0);
            progressCallback(currentIteration, maxK,
                "Calculating WPI and selecting prevalent patterns...",
                progressPercent);
        }

        // Step 11: calculate_WPI(Ck, Tk)
        // Step 12: Pk = select_prevalent_patterns(Ck, Tk, min_prev)
        auto t4_start = std::chrono::high_resolution_clock::now();
        prevColocations = selectPrevColocations(
            fiteredCandidates,
            tableInstances,
            minPrev
        );
        auto t4_end = std::chrono::high_resolution_clock::now();
        printDuration("Step 11-12: select_prevalent_patterns (k=" + std::to_string(k) + ")", t4_start, t4_end);


        if (!prevColocations.empty()) {
            allPrevalentColocations.insert(
                allPrevalentColocations.end(),
                prevColocations.begin(),
                prevColocations.end()
            );

            std::cout << "[DEBUG] Step 12: Found " << prevColocations.size() << " prevalent patterns for k=" << k << "\n";

            if (progressCallback) {
                double progressPercent = std::min(95.0, (static_cast<double>(currentIteration) / maxK) * 95.0);
                progressCallback(currentIteration, maxK,
                    "Found " + std::to_string(prevColocations.size()) + " prevalent k=" + std::to_string(k) + " co-locations",
                    progressPercent);
            }
        }
        else {
            if (progressCallback) {
                double progressPercent = std::min(95.0, (static_cast<double>(currentIteration) / maxK) * 95.0);
                progressCallback(currentIteration, maxK,
                    "No prevalent k=" + std::to_string(k) + " co-locations found",
                    progressPercent);
            }
        }

        // Prepare Tk-1 for next iteration
        prevTableInstances = std::move(tableInstances);

        // Step 13: k = k + 1
        k++;
    } // Step 14: end while

    // Step 15: return union(P2, ..., Pk-1)
    if (progressCallback) {
        progressCallback(maxK, maxK,
            "Mining completed! Total prevalent co-locations: " + std::to_string(allPrevalentColocations.size()),
            100.0);
    }

    auto minerEnd = std::chrono::high_resolution_clock::now();
    printDuration("TOTAL MINING TIME (Algorithm 1)", minerStart, minerEnd);

    return allPrevalentColocations;
}



std::vector<Colocation> JoinlessMiner::generateCandidates(
    const std::vector<Colocation>& prevPrevalent,
    const std::map<FeatureType, int>& featureCount)
{
    // Implementation of gen_candidate_patterns (Step 8)
    std::vector<Colocation> candidates;

    if (prevPrevalent.empty()) {
        return candidates;
    }

    size_t patternSize = prevPrevalent[0].size();

    std::set<Colocation> prevSet(prevPrevalent.begin(), prevPrevalent.end());

    // Join phase: generate k-size candidates from (k-1)-size prevalent patterns
    for (size_t i = 0; i < prevPrevalent.size(); i++) {
        for (size_t j = i + 1; j < prevPrevalent.size(); j++) {
            // Take prefix of k-1 first element
            Colocation prefix1(prevPrevalent[i].begin(),
                prevPrevalent[i].end() - 1);
            Colocation prefix2(prevPrevalent[j].begin(),
                prevPrevalent[j].end() - 1);

            // Just join when the prefix is equal
            if (prefix1 != prefix2) {
                continue;
            }
            
            // Generate new candidate
			std::set<FeatureType> candidateSet;
            if (featureCount.at(prevPrevalent[i].back()) <= featureCount.at(prevPrevalent[j].back())) {
                candidateSet = std::set<FeatureType>(prevPrevalent[i].begin(),
                    prevPrevalent[i].end());
                candidateSet.insert(prevPrevalent[j].back());
            }else {
                candidateSet = std::set<FeatureType>(prevPrevalent[j].begin(),
                    prevPrevalent[j].end());
                candidateSet.insert(prevPrevalent[i].back());
            }
            
            if (candidateSet.size() != patternSize + 1) {
                continue;
            }

            Colocation candidate(candidateSet.begin(), candidateSet.end());
        }
    }

    // Remove duplication
    std::sort(candidates.begin(), candidates.end());
    candidates.erase(std::unique(candidates.begin(), candidates.end()),
        candidates.end());

    return candidates;
}


std::vector<Colocation> JoinlessMiner::filterCandidates(
    const std::vector<Colocation>& candidates,
    const std::vector<Colocation>& prevPrevalent,
	const std::vector<ColocationInstance>& tableInstance,
    double minPrev,
    std::map<FeatureType, int> featureCount,
    double delta)
{
    // Implementation of filter_candidate_patterns (Step 9)
    std::vector<Colocation> filteredCandidates;
    if (candidates.empty() || prevPrevalent.empty()) {
        return filteredCandidates;
    }
    std::set<Colocation> prevSet(prevPrevalent.begin(), prevPrevalent.end());
    for (const auto& candidate : candidates) {
        bool isValid = true;
        // Generate all (k-1)-size subsets
        for (size_t i = 0; i < candidate.size(); i++) {
            Colocation subset;
            for (size_t j = 0; j < candidate.size(); j++) {
                if (j != i) subset.push_back(candidate[j]);
            }

            // --- CASE 1: Subset contains f_min (Lemma 2) ---
            // If we removed an element at index i != 0, the subset still keeps candidate[0] (f_min).
            if (i != 0) {
                // Lemma 2: If a subset containing f_min is NOT prevalent, C is not prevalent 
                if (prevSet.find(subset) == prevSet.end()) {
                    isValid = false;
                    break; // Prune immediately
                }
            }

            // --- CASE 2: Subset does NOT contain f_min (Lemma 3) ---
            // This happens when i == 0 (we removed f_min). Subset is {f2, ..., fk}.
            else {
                // Lemma 3: Check upper bound condition 
                // Condition: PI(subset) * w(f_max, C) < min_prev => Prune

                // 1. Find f_max (feature with max instances in C)
                // Assuming sorted input, f_max is the last element
                FeatureType f_max = candidate.back();

                // 2. Calculate Weight w(f_max, C) = 1 / RI(f_max, C) 
                double RI = calculateRareIntensity(f_max, candidate, featureCount, delta);
                double w = 1.0 / RI;

                // 3. Get PI of the subset (needs to be looked up from previous results)
                double piSubset = calculatePI(subset, tableInstance, featureCount);

                // Check Lemma 3 inequality
                if (piSubset * w < minPrev) {
                    isValid = false;
                    break; // Prune
                }
            }
        }
        if (isValid) {
            filteredCandidates.push_back(candidate);
        }
    }
    return filteredCandidates;
}