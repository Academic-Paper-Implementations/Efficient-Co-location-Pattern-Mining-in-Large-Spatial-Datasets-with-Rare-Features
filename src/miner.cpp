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
    NRTree orderedNRTree,
    const std::vector<SpatialInstance>& instances,
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
    std::map<FeatureType, int> featureCount = countInstancesByFeature(instances);

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
        std::vector<Colocation> candidates = generateCandidates(prevColocations);
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
        std::vector<Colocation> fiteredCandidates = filterCandidates(candidates, prevColocations);
        auto t2_end = std::chrono::high_resolution_clock::now();
        printDuration("Step 9: filter_candidate_patterns (k=" + std::to_string(k) + ")", t2_start, t2_end);

        // Step 10: Tk = gen_table_instances(Ck, Tk-1, ordered-NR-tree)
        auto t3_start = std::chrono::high_resolution_clock::now();
        std::vector<ColocationInstance> tableInstances = genTableInstance(
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
    const std::vector<Colocation>& prevPrevalent)
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
            std::set<FeatureType> candidateSet(prevPrevalent[i].begin(),
                prevPrevalent[i].end());
            candidateSet.insert(prevPrevalent[j].back());

            if (candidateSet.size() != patternSize + 1) {
                continue;
            }

            Colocation candidate(candidateSet.begin(), candidateSet.end());

            // Prune phase: check subsets (Lemma 2)
            bool allSubsetsValid = true;
            std::vector<FeatureType> candFeatures = candidate;

            for (size_t idx = 0; idx < candFeatures.size(); idx++) {
                Colocation subset = candFeatures;
                subset.erase(subset.begin() + idx);

                if (prevSet.find(subset) == prevSet.end()) {
                    allSubsetsValid = false;
                    break;
                }
            }

            if (allSubsetsValid) {
                candidates.push_back(candidate);
            }
        }
    }

    // Remove duplication
    std::sort(candidates.begin(), candidates.end());
    candidates.erase(std::unique(candidates.begin(), candidates.end()),
        candidates.end());

    return candidates;
}


std::vector<ColocationInstance> JoinlessMiner::filterStarInstances(
    const std::vector<Colocation>& candidates,
    const std::pair<FeatureType, std::vector<StarNeighborhood>>& starNeigh)
{
    // Part of row instance generation logic
    std::vector<ColocationInstance> filteredInstances;
    FeatureType centerType = starNeigh.first;

    // Filter candidates to only those with this center type as first element
    std::vector<const Colocation*> relevantCandidates;
    for (const auto& cand : candidates) {
        if (!cand.empty() && cand[0] == centerType) {
            relevantCandidates.push_back(&cand);
        }
    }

    if (relevantCandidates.empty()) return filteredInstances;

    // Iterate through each star neighborhood (concept from Ordered NR-Tree branch)
    for (const auto& star : starNeigh.second) {

        // Build a map of neighbors by feature type for fast lookup
        // Using const pointer since star is const
        std::unordered_map<FeatureType, std::vector<const SpatialInstance*>> neighborMap;

        for (auto neighbor : star.neighbors) {
            neighborMap[neighbor->type].push_back(neighbor);
        }

        // Check each relevant candidate pattern
        for (const auto* candPtr : relevantCandidates) {
            const auto& candidate = *candPtr;

            std::vector<const SpatialInstance*> currentInstance;
            currentInstance.reserve(candidate.size());

            // Add center instance as first element
            currentInstance.push_back(star.center);

            // Recursive function to find combinations (row instances)
            findCombinations(candidate, 1, currentInstance, neighborMap, filteredInstances);
        }
    }

    return filteredInstances;
}


std::vector<ColocationInstance> JoinlessMiner::filterCliqueInstances(
    const std::vector<Colocation>& candidates,
    const std::vector<ColocationInstance>& instances,
    const std::vector<ColocationInstance>& prevInstances
) {
    // Helper logic for Step 10 (gen_table_instances) validation
    // ========================================================================
    // PREPARE LOOKUP STRUCTURES
    // ========================================================================

    // 1. Create a set of valid candidate patterns for quick lookup
    std::set<Colocation> validCandidatePatterns(candidates.begin(), candidates.end());

    // 2. Create a set of previous instances for quick lookup (Check against Tk-1)
    std::set<std::vector<std::string>> validPrevIds;
    for (const auto& prevInst : prevInstances) {
        std::vector<std::string> ids;
        ids.reserve(prevInst.size());
        for (const auto* ptr : prevInst) {
            ids.push_back(ptr->id);
        }
        validPrevIds.insert(ids);
    }

    // ========================================================================
    // PREPARE THREAD BUFFERS
    // ========================================================================
    int num_threads = omp_get_max_threads();
    std::vector<std::vector<ColocationInstance>> thread_buffers(num_threads);

    // ========================================================================
    // PARALLEL FILTERING
    // ========================================================================
#pragma omp parallel for
    for (size_t i = 0; i < instances.size(); ++i) {
        const auto& instance = instances[i];
        int thread_id = omp_get_thread_num();

        // Safety check
        if (instance.size() < 2) continue;

        Colocation currentPattern;
        currentPattern.reserve(instance.size());
        for (const auto* ptr : instance) {
            currentPattern.push_back(ptr->type);
        }

        // If current pattern is not a valid candidate, skip
        if (validCandidatePatterns.find(currentPattern) == validCandidatePatterns.end()) {
            continue;
        }

        std::vector<std::string> subInstanceIds;
        subInstanceIds.reserve(instance.size() - 1);

        // Generate (k-1)-subset by removing the first instance
        for (size_t j = 1; j < instance.size(); ++j) {
            subInstanceIds.push_back(instance[j]->id);
        }

        // Check if the (k-1)-subset exists in previous instances (Tk-1)
        if (validPrevIds.find(subInstanceIds) != validPrevIds.end()) {
            thread_buffers[thread_id].push_back(instance);
        }
    }

    // ========================================================================
    // COMBINE RESULTS
    // ========================================================================
    std::vector<ColocationInstance> filteredInstances;
    for (const auto& buffer : thread_buffers) {
        filteredInstances.insert(filteredInstances.end(), buffer.begin(), buffer.end());
    }

    return filteredInstances;
}


std::vector<Colocation> JoinlessMiner::selectPrevColocations(
    const std::vector<Colocation>& candidates,
    const std::vector<ColocationInstance>& instances,
    double minPrev,
    const std::map<FeatureType, int>& featureCount)
{
    // Implementation of Step 11: calculate_WPI and Step 12: select_prevalent_patterns
    std::vector<Colocation> coarsePrevalent;

    // ========================================================================
    // Data structure for aggregation (calculating participation)
    // ========================================================================
    // Key: Candidate (Pattern)
    // Value: Map<FeatureType, Set<InstanceID>> - count unique instances for each feature
    std::map<Colocation, std::map<FeatureType, std::set<std::string>>> candidateStats;

    // Initialize stats map for all candidates
    for (const auto& cand : candidates) {
        candidateStats[cand]; // Create empty entry
    }

    // ========================================================================
    // Single pass through instances (Tk)
    // ========================================================================
    for (const ColocationInstance& instance : instances) {
        // Extract pattern from instance
        Colocation patternKey;
        for (const auto& instPtr : instance) {
            patternKey.push_back(instPtr->type);
        }

        // Check if this pattern is in the candidates of interest
        auto it = candidateStats.find(patternKey);
        if (it != candidateStats.end()) {
            // Update participating instances statistics
            for (const auto& instPtr : instance) {
                it->second[instPtr->type].insert(instPtr->id);
            }
        }
    }

    // ========================================================================
    // Calculate Metric and filter
    // ========================================================================
    for (const auto& item : candidateStats) {
        const Colocation& candidate = item.first;
        const auto& participatingMap = item.second; // Map<Feature, Set<ID>>

        double min_participation_ratio = 1.0;
        bool possible = true;

        // Loop through features to calculate PR / WPI components
        for (const auto& feature : candidate) {
            // Get total global instance count for this feature
            auto totalIt = featureCount.find(feature);
            if (totalIt == featureCount.end() || totalIt->second == 0) {
                possible = false;
                break;
            }
            double totalFeatureCount = (double)(totalIt->second);

            // Get count of instances participating in pattern
            int participatedCount = 0;
            auto partIt = participatingMap.find(feature);
            if (partIt != participatingMap.end()) {
                participatedCount = partIt->second.size();
            }

            // [ANOMALY FLAG]: Code calculates standard Participation Ratio (PR).
            // Paper requires Weighted Participation Ratio (WPR) = PR * w(f, C).
            double ratio = (double)participatedCount / totalFeatureCount;
            if (ratio < min_participation_ratio) {
                min_participation_ratio = ratio;
            }
        }

        // Check against min_prev
        // Note: For WPI, this should clearly verify WPI(C) >= min_prev
        if (possible && min_participation_ratio >= minPrev) {
            coarsePrevalent.push_back(candidate);
        }
    }

    return coarsePrevalent;
}