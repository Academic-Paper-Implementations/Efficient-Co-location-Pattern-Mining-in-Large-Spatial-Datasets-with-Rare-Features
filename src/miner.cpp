/**
 * @file miner.cpp
 * @brief Implementation of the ordered NR-tree based mining algorithm
 * * This file implements the core mining algorithm that discovers prevalent spatial
 * co-location patterns using the Ordered NR-Tree structure as described in Algorithm 1.
 */

#include "miner.h"
#include "utils.h"
#include "constants.h"
#include "neighborhood_mgr.h"
#include "NRTree.h"
#include "types.h"
#include <algorithm>
#include <set>
#include <map>
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
    std::map<Colocation, std::vector<ColocationInstance>> prevTableInstances;
    for (const auto& instance : instances) {
        Colocation key = { instance.type };
        ColocationInstance row = { &instance };
        prevTableInstances[key].push_back(row);
    }

    std::vector<Colocation> allPrevalentColocations;

    // Estimate total iterations (max pattern size is number of types)
    int maxK = static_cast<int>(types.size());
    int currentIteration = 0;
    int totalIterations = 0;

    if (progressCallback) {
        progressCallback(0, maxK, "Initializing mining process (Steps 1-6)...", 0.0);
    }

    // Step 6: let P1 = F, T1 = O, k = 2
    for (auto t : sortedTypes) prevColocations.push_back({ t });

    // Step 7: while Pk-1 not empty do
    while (!prevColocations.empty()) {
        currentIteration++;
        totalIterations = currentIteration;

        // Calculate progress
        double progressPercent = std::min(Constants::MAX_PROGRESS_PERCENT, (static_cast<double>(currentIteration) / maxK) * Constants::MAX_PROGRESS_PERCENT);

        if (progressCallback) {
            progressCallback(currentIteration, maxK,
                "Processing iteration k=" + std::to_string(k) + "...",
                progressPercent);
        }
        std::map<Colocation, std::vector<ColocationInstance>> tableInstances;

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
            double progressPercent = std::min(Constants::MAX_PROGRESS_PERCENT, (static_cast<double>(currentIteration) / maxK) * Constants::MAX_PROGRESS_PERCENT);
            progressCallback(currentIteration, maxK,
                "Filtering candidates (Lemma 2 & 3)...",
                progressPercent);
        }

        // Step 9: filter_candidate_patterns(Ck, Pk-1)
        // Uses Lemma 2 and Lemma 3 to prune search space
        auto t2_start = std::chrono::high_resolution_clock::now();
		std::vector<Colocation> fiteredCandidates = candidates;
        if (k!=2){
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
            double progressPercent = std::min(Constants::MAX_PROGRESS_PERCENT, (static_cast<double>(currentIteration) / maxK) * Constants::MAX_PROGRESS_PERCENT);
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
            minPrev,
			featureCount,
			delta
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
                double progressPercent = std::min(Constants::MAX_PROGRESS_PERCENT, (static_cast<double>(currentIteration) / maxK) * Constants::MAX_PROGRESS_PERCENT);
                progressCallback(currentIteration, maxK,
                    "Found " + std::to_string(prevColocations.size()) + " prevalent k=" + std::to_string(k) + " co-locations",
                    progressPercent);
            }
        }
        else {
            if (progressCallback) {
                double progressPercent = std::min(Constants::MAX_PROGRESS_PERCENT, (static_cast<double>(currentIteration) / maxK) * Constants::MAX_PROGRESS_PERCENT);
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
			Colocation candidate;
            if (featureCount.at(prevPrevalent[i].back()) <= featureCount.at(prevPrevalent[j].back())) {
                candidate = prevPrevalent[i];
                candidate.push_back(prevPrevalent[j].back());
            }else {
				candidate = prevPrevalent[j];
				candidate.push_back(prevPrevalent[i].back());
            }
            
            candidates.push_back(candidate);
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
	const std::map<Colocation, std::vector<ColocationInstance>>& tableInstance,
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

// Helper function to find neighbors of an instance for a specific feature type from NRTree
// Returns Neigh(o, f) - all neighbors of instance o that have feature type f
std::vector<const SpatialInstance*> JoinlessMiner::findNeighbors(
    const NRTree& tree,
    const SpatialInstance* instance,
    const FeatureType& featureType
) {
    std::vector<const SpatialInstance*> result;
    const NRNode* root = tree.getRoot();
    if (!root) return result;

    // Traverse NRTree structure:
    // Level 1: FEATURE_NODE (center features)
    // Level 2: INSTANCE_NODE (center instances)
    // Level 3: FEATURE_NODE (neighbor features)
    // Level 4: INSTANCE_VECTOR_NODE (neighbor instances)

    // Find the feature node for the instance's feature type (Level 1)
    for (const auto* featureNode : root->children) {
        if (featureNode->type == FEATURE_NODE && featureNode->featureType == instance->type) {
            // Find the instance node for this specific instance (Level 2)
            for (const auto* instanceNode : featureNode->children) {
                if (instanceNode->type == INSTANCE_NODE && instanceNode->data == instance) {
                    // Find the neighbor feature node for the target feature type (Level 3)
                    for (const auto* neighborFeatureNode : instanceNode->children) {
                        if (neighborFeatureNode->type == FEATURE_NODE && 
                            neighborFeatureNode->featureType == featureType) {
                            // Get the instance vector node (Level 4)
                            if (!neighborFeatureNode->children.empty()) {
                                const auto* instanceVectorNode = neighborFeatureNode->children[0];
                                if (instanceVectorNode->type == INSTANCE_VECTOR_NODE) {
                                    // Return all neighbor instances
                                    return instanceVectorNode->instanceVector;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    return result;
}

// Helper function to calculate S(I, f) = Neigh(o1, f) ∩ ··· ∩ Neigh(ok, f) (Definition 8)
// Returns the intersection of neighbors for all instances in I with feature type f
std::vector<const SpatialInstance*> JoinlessMiner::findExtendedSet(
    const NRTree& tree,
    const ColocationInstance& instance,
    const FeatureType& featureType
) {
    if (instance.empty()) {
        return std::vector<const SpatialInstance*>();
    }

    // Start with neighbors of the first instance
    std::vector<const SpatialInstance*> intersection = findNeighbors(tree, instance[0], featureType);
    
    // Intersect with neighbors of remaining instances
    for (size_t i = 1; i < instance.size(); i++) {
        std::vector<const SpatialInstance*> neighbors = findNeighbors(tree, instance[i], featureType);
        
        // Calculate intersection: keep only instances that are in both sets
        std::vector<const SpatialInstance*> newIntersection;
        std::set<const SpatialInstance*> neighborSet(neighbors.begin(), neighbors.end());
        
        for (const auto* inst : intersection) {
            if (neighborSet.find(inst) != neighborSet.end()) {
                newIntersection.push_back(inst);
            }
        }
        
        intersection = std::move(newIntersection);
        
        // Early termination if intersection becomes empty
        if (intersection.empty()) {
            break;
        }
    }

    return intersection;
}

//// Main function to generate table instances for candidate patterns (Step 10)
//// Uses Definition 8 and Lemma 4 to extend (k-1)-size instances to k-size instances
//std::map<Colocation, std::vector<ColocationInstance>> JoinlessMiner::genTableInstance(
//    const std::vector<Colocation>& candidates,
//    const std::map<Colocation, std::vector<ColocationInstance>>& prevTableInstances,
//    const NRTree& orderedNRTree
//) {
//    std::map<Colocation, std::vector<ColocationInstance>> result;
//
//    // For each candidate pattern C of size k
//    for (const auto& candidate : candidates) {
//        std::vector<ColocationInstance> tableInstances;
//
//        // For each table instance I of size k-1
//        for (const auto& prevInstance : prevTableInstances) {
//            // Check if prevInstance can be extended to form an instance of candidate
//            // According to the algorithm, candidates are generated from (k-1)-size patterns
//            // using Apriori-gen, so prevInstance should match a (k-1)-size subset of candidate
//            
//            // Check if prevInstance size matches (k-1)
//            if (prevInstance.size() != candidate.size() - 1) {
//                continue;
//            }
//
//            // Build set of features in prevInstance
//            std::set<FeatureType> prevInstanceFeatures;
//            for (const auto* inst : prevInstance) {
//                prevInstanceFeatures.insert(inst->type);
//            }
//
//            // Build set of features in candidate
//            std::set<FeatureType> candidateFeatures(candidate.begin(), candidate.end());
//
//            // Check if prevInstance is a subset of candidate (exactly k-1 features match)
//            if (prevInstanceFeatures.size() != candidateFeatures.size() - 1) {
//                continue;
//            }
//
//            // Find the new feature f (the one in candidate but not in prevInstance)
//            FeatureType newFeature = "";
//            bool isSubset = true;
//            
//            for (const auto& feature : candidateFeatures) {
//                if (prevInstanceFeatures.find(feature) == prevInstanceFeatures.end()) {
//                    if (newFeature.empty()) {
//                        newFeature = feature;
//                    } else {
//                        // More than one feature missing - prevInstance doesn't match
//                        isSubset = false;
//                        break;
//                    }
//                }
//            }
//
//            // Verify that all features in prevInstance are in candidate
//            for (const auto& feature : prevInstanceFeatures) {
//                if (candidateFeatures.find(feature) == candidateFeatures.end()) {
//                    isSubset = false;
//                    break;
//                }
//            }
//
//            if (!isSubset || newFeature.empty()) {
//                continue;
//            }
//
//            // Now we need to reorder prevInstance to match the order in candidate
//            // Create a mapping from feature to instance in prevInstance
//            std::map<FeatureType, const SpatialInstance*> featureToInstance;
//            for (const auto* inst : prevInstance) {
//                featureToInstance[inst->type] = inst;
//            }
//
//            // Build ordered instance matching candidate order (excluding newFeature)
//            ColocationInstance orderedPrevInstance;
//            for (const auto& feature : candidate) {
//                if (feature != newFeature) {
//                    auto it = featureToInstance.find(feature);
//                    if (it != featureToInstance.end()) {
//                        orderedPrevInstance.push_back(it->second);
//                    } else {
//                        // This shouldn't happen if isSubset is true, but check anyway
//                        isSubset = false;
//                        break;
//                    }
//                }
//            }
//
//            if (!isSubset || orderedPrevInstance.size() != candidate.size() - 1) {
//                continue;
//            }
//
//            // Calculate S(I, f) = Neigh(o1, f) ∩ ··· ∩ Neigh(ok, f) (Definition 8)
//            std::vector<const SpatialInstance*> extendedSet = findExtendedSet(
//                orderedNRTree, orderedPrevInstance, newFeature
//            );
//
//            // According to Lemma 4: if S(I, f) ≠ Ø, append any instance o from S(I, f) to I
//            // to form a row instance I' = {o1, ..., ok, o} of co-location C' = {f1, ..., fk, f}
//            for (const auto* newInstance : extendedSet) {
//                // Verify that the new instance has the correct feature type
//                if (newInstance->type != newFeature) {
//                    continue;
//                }
//
//                // Create new row instance I' = orderedPrevInstance ∪ {o}
//                // The order should match candidate: [f1-instance, ..., fk-instance]
//                ColocationInstance newRowInstance = orderedPrevInstance;
//                newRowInstance.push_back(newInstance);
//                tableInstances.push_back(newRowInstance);
//            }
//        }
//
//        // Store table instances for this candidate pattern
//        if (!tableInstances.empty()) {
//            result[candidate] = tableInstances;
//        }
//    }
//
//    return result;
//}


std::map<Colocation, std::vector<ColocationInstance>> JoinlessMiner::genTableInstance(
    const std::vector<Colocation>& candidates,
    const std::map<Colocation, std::vector<ColocationInstance>>& prevTableInstances,
    const NRTree& orderedNRTree
) {
    std::map<Colocation, std::vector<ColocationInstance>> result;

    // Iterate through each candidate pattern C of size k
    for (const auto& candidate : candidates) {
        if (candidate.empty()) continue;

        // 1. Split candidate into prefix (k-1 features) and the new feature
        // Assumes candidate features are sorted, so prefix is the first k-1 elements
        Colocation subPattern(candidate.begin(), candidate.end() - 1);
        FeatureType newFeature = candidate.back();

        // 2. Fast lookup: Find existing instances of the prefix pattern
        auto it = prevTableInstances.find(subPattern);
        if (it == prevTableInstances.end()) {
            continue; // No instances to extend
        }

        std::vector<ColocationInstance> newTableRows;
        // Access the vector of instances directly (Value of the Map)
        const std::vector<ColocationInstance>& prevInstancesList = it->second;

        // 3. Try to extend each existing instance I with the new feature f
        for (const auto& prevInstance : prevInstancesList) {

            // Calculate intersection of neighbors using NRTree: S(I, f)
            std::vector<const SpatialInstance*> extendedSet = findExtendedSet(
                orderedNRTree, prevInstance, newFeature
            );

            // 4. Create new instances I' = I + {o}
            for (const auto* neighbor : extendedSet) {
                // strict type check
                if (neighbor->type == newFeature) {
                    ColocationInstance newRow = prevInstance;
                    newRow.push_back(neighbor);
                    newTableRows.push_back(newRow);
                }
            }
        }

        // Store results if any instances were found
        if (!newTableRows.empty()) {
            result[candidate] = newTableRows;
        }
    }

    return result;
}


std::vector<Colocation> JoinlessMiner::selectPrevColocations(
    const std::vector<Colocation>& candidates,
    const std::map<Colocation, std::vector<ColocationInstance>>& tableInstances,
    double minPrev,
    const std::map<FeatureType, int>& featureCount,
    double delta
) {
    std::vector<Colocation> prevalentPatterns;

    for (const auto& candidate : candidates) {
        // Init WPI with a safe max value
        double wpi = 1.0;
        bool first = true;

        for (const FeatureType& feature : candidate) {
            double pr = calculatePR(feature, candidate, tableInstances, featureCount);
            double ri = calculateRareIntensity(feature, candidate, featureCount, delta);

            // Safety check for RI
            double w = 0.0;
            if (ri > Constants::EPSILON_SMALL) { // Epsilon check
                w = 1.0 / ri;
            }
            else {
                // If RI is 0 (should imply feature not in pattern or error), weight is huge or handled
                w = 0.0;
            }

            double wpr = pr * w;

            if (first) {
                wpi = wpr;
                first = false;
            }
            else {
                if (wpr < wpi) {
                    wpi = wpr;
                }
            }
        }

        // Check threshold
        if (wpi >= minPrev) {
            prevalentPatterns.push_back(candidate);
        }
    }

    return prevalentPatterns;
}
