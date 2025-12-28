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
#include <unordered_set>
#include <map>
#include <string>
#include <iostream>
#include <omp.h> 
#include <iomanip>
#include <chrono>
#include <unordered_set>

std::vector<Colocation> JoinlessMiner::mineColocations(
    double minPrev,
    NRTree& orderedNRTree,
    const std::vector<SpatialInstance>& instances,
    const std::map<FeatureType, int>& featureCount,
    ProgressCallback progressCb
) {
    auto minerStart = std::chrono::high_resolution_clock::now();
    this->progressCallback = progressCb;

    // --- INIT ---
    int k = 2;
    std::vector<FeatureType> types = getAllObjectTypes(instances);
    std::vector<FeatureType> sortedTypes = featureSort(types, instances);
    double delta = calculateDelta(sortedTypes, featureCount);

    // DEBUG: In ra thứ tự Feature sau khi sort (Rất quan trọng)
    std::cout << "\n[DEBUG INIT] Feature Sort Order (Rare -> Frequent):\n";
    for (const auto& t : sortedTypes) {
        std::cout << "   " << t << ": " << featureCount.at(t) << " instances\n";
    }
    std::cout << "------------------------------------------------\n";

    std::vector<Colocation> prevColocations;
    std::map<Colocation, std::vector<ColocationInstance>> prevTableInstances;

    // Khởi tạo T1 (Table Instance k=1)
    // Map: {Key: [FeatureType], Value: List of rows}
    for (const auto& instance : instances) {
        Colocation key = { instance.type };
        ColocationInstance row = { &instance };
        prevTableInstances[key].push_back(row);
    }

    std::vector<Colocation> allPrevalentColocations;

    // Khởi tạo P1 (Prevalent k=1)
    for (auto t : sortedTypes) prevColocations.push_back({ t });

    // --- MAIN LOOP ---
    while (!prevColocations.empty()) {
        std::cout << "\n>>> ITERATION K = " << k << " <<<\n";

        std::map<Colocation, std::vector<ColocationInstance>> tableInstances;

        // 1. Generate Candidates
        auto t1_start = std::chrono::high_resolution_clock::now();
        std::vector<Colocation> candidates = generateCandidates(prevColocations, featureCount);
        auto t1_end = std::chrono::high_resolution_clock::now();

        std::cout << "   [GEN] Generated: " << candidates.size() << " candidates. ";
        printDuration("", t1_start, t1_end);

        // --- DEBUG CHI TIẾT CANDIDATES ---
        if (!candidates.empty()) {
            std::cout << "      List of Candidates:\n";
            for (const auto& cand : candidates) {
                std::cout << "      - { ";
                for (size_t i = 0; i < cand.size(); ++i) std::cout << cand[i] << (i < cand.size() - 1 ? ", " : "");
                std::cout << " }\n";
            }
        }
        // --------------------------------

        if (candidates.empty()) break;

        // 2. Filter Candidates
        auto t2_start = std::chrono::high_resolution_clock::now();
        std::vector<Colocation> filteredCandidates = candidates;
        if (k != 2) {
            filteredCandidates = filterCandidates(candidates, prevColocations, prevTableInstances, minPrev, featureCount, delta);
        }
        auto t2_end = std::chrono::high_resolution_clock::now();

        std::cout << "   [FLT] Filtered: " << candidates.size() << " -> " << filteredCandidates.size() << " candidates. ";
        printDuration("", t2_start, t2_end);

        if (filteredCandidates.empty()) break;

        // 3. Generate Table Instances (Phần quan trọng nhất cần check)
        auto t3_start = std::chrono::high_resolution_clock::now();
        tableInstances = genTableInstance(filteredCandidates, prevTableInstances, orderedNRTree);
        auto t3_end = std::chrono::high_resolution_clock::now();

        size_t totalRows = 0;
        for (const auto& pair : tableInstances) totalRows += pair.second.size();

        std::cout << "   [VER] Verified: Built instance table with " << totalRows << " rows. ";
        printDuration("", t3_start, t3_end);

        // --- DEBUG CHI TIẾT INSTANCE TABLE ---
        std::cout << "      [DEBUG VERIFY] Instance Table Details:\n";
        if (tableInstances.empty()) {
            std::cout << "      !!! WARNING: No instances found for ANY candidate.\n";
            std::cout << "      !!! Check: 1. NRTree correctness. 2. NeighborDistance (d) too small?\n";
        }
        else {
            for (const auto& cand : filteredCandidates) {
                // In ra candidate
                std::cout << "      Pattern { ";
                for (const auto& f : cand) std::cout << f << " ";
                std::cout << "}: ";

                // Kiểm tra xem candidate này có instance nào không
                auto it = tableInstances.find(cand);
                if (it != tableInstances.end()) {
                    const auto& rows = it->second;
                    std::cout << rows.size() << " instances found.\n";

                    // In thử tối đa 3 dòng đầu tiên để check ID
                    int printLimit = 3;
                    for (size_t i = 0; i < std::min((size_t)rows.size(), (size_t)printLimit); ++i) {
                        std::cout << "         Row " << i + 1 << ": [ ";
                        for (const auto* inst : rows[i]) {
                            std::cout << inst->id << "(" << inst->type << ") ";
                        }
                        std::cout << "]\n";
                    }
                    if (rows.size() > printLimit) std::cout << "         ... (more)\n";

                }
                else {
                    std::cout << "0 instances (Empty).\n";
                }
            }
        }
        // -------------------------------------

        // 4. Select Prevalent
        auto t4_start = std::chrono::high_resolution_clock::now();
        prevColocations = selectPrevColocations(
            filteredCandidates,
            tableInstances,
            minPrev,
            featureCount,
            delta
        );
        auto t4_end = std::chrono::high_resolution_clock::now();

        std::cout << "   [SEL] Selected: " << prevColocations.size() << " prevalent patterns. ";
        printDuration("", t4_start, t4_end);

        // Debug PI values (nếu cần thiết thì bỏ comment đoạn này)
        /*
        for (const auto& cand : filteredCandidates) {
             double pi = calculatePI(cand, tableInstances, featureCount);
             std::cout << "      PI({ ";
             for(auto f : cand) std::cout << f << " ";
             std::cout << "}) = " << pi << " (Threshold: " << minPrev << ")\n";
        }
        */

        if (!prevColocations.empty()) {
            allPrevalentColocations.insert(allPrevalentColocations.end(), prevColocations.begin(), prevColocations.end());
        }

        prevTableInstances = std::move(tableInstances);
        k++;
    }

    auto minerEnd = std::chrono::high_resolution_clock::now();
    std::cout << "\n[MINER] Total Mining Time: "
        << std::chrono::duration<double, std::milli>(minerEnd - minerStart).count() << " ms\n";

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
    std::vector<Colocation> sortedPrev = prevPrevalent;
    std::sort(sortedPrev.begin(), sortedPrev.end());
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
                if (!std::binary_search(sortedPrev.begin(), sortedPrev.end(), subset)) {
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
                if (instanceNode->type == INSTANCE_NODE && instanceNode->data->id == instance->id) {
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
        return {};
    }

    // Start with neighbors of the first instance
    std::vector<const SpatialInstance*> intersection = findNeighbors(tree, instance[0], featureType);

    // Intersect with neighbors of remaining instances
    for (size_t i = 1; i < instance.size(); i++) {
        std::vector<const SpatialInstance*> neighbors = findNeighbors(tree, instance[i], featureType);

        if (neighbors.empty()) {
            return {}; // Một ông không có neighbor thì giao bằng rỗng luôn
        }

        // --- [SỬA TẠI ĐÂY] ---
        // Dùng ID để tạo bảng tra cứu (Lookup Table)
        // Giả sử ID là int hoặc string đều dùng được cách này
        std::set<decltype(instance[0]->id)> neighborIDs;
        for (const auto* ptr : neighbors) {
            neighborIDs.insert(ptr->id);
        }

        std::vector<const SpatialInstance*> newIntersection;

        // Chỉ giữ lại những thằng trong intersection có ID nằm trong neighborIDs
        for (const auto* ptr : intersection) {
            if (neighborIDs.count(ptr->id)) {
                newIntersection.push_back(ptr);
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


// Hàm helper để in pattern (VD: {A, B})
void printDebugPattern(const Colocation& col) {
    std::cout << "{";
    for (size_t i = 0; i < col.size(); ++i) {
        std::cout << col[i] << (i < col.size() - 1 ? ", " : "");
    }
    std::cout << "}";
}

std::map<Colocation, std::vector<ColocationInstance>> JoinlessMiner::genTableInstance(
    const std::vector<Colocation>& candidates,
    const std::map<Colocation, std::vector<ColocationInstance>>& prevTableInstances,
    const NRTree& orderedNRTree
) {
    std::map<Colocation, std::vector<ColocationInstance>> result;

    // Iterate through each candidate pattern C of size k
    for (const auto& candidate : candidates) {

        // --- [DEBUG 1] Kiểm tra candidate rỗng ---
        if (candidate.empty()) {
            std::cout << "[DEBUG] SKIP: Candidate is empty.\n";
            continue;
        }

        // 1. Split candidate into prefix (k-1 features) and the new feature
        Colocation subPattern(candidate.begin(), candidate.end() - 1);
        FeatureType newFeature = candidate.back();

        // 2. Fast lookup: Find existing instances of the prefix pattern
        auto it = prevTableInstances.find(subPattern);

        // --- [DEBUG 2] Kiểm tra Prefix (quan trọng nhất cho lỗi size 2) ---
        if (it == prevTableInstances.end()) {
            std::cout << "[DEBUG] SKIP Candidate ";
            printDebugPattern(candidate);
            std::cout << ". Reason: Prefix ";
            printDebugPattern(subPattern);
            std::cout << " NOT FOUND in prevTableInstances.\n";

            // Gợi ý lỗi cụ thể nếu đang chạy k=2
            if (candidate.size() == 2) {
                std::cout << "    -> HINT: For k=2, prevTableInstances must contain size-1 patterns (e.g., {'A'}). Did you initialize T1 correctly?\n";
            }
            continue;
        }

        std::vector<ColocationInstance> newTableRows;
        const std::vector<ColocationInstance>& prevInstancesList = it->second;

        // --- [DEBUG 3] Prefix tồn tại nhưng không có instance nào ---
        if (prevInstancesList.empty()) {
            std::cout << "[DEBUG] SKIP Candidate ";
            printDebugPattern(candidate);
            std::cout << ". Reason: Prefix found but has 0 instances.\n";
            continue;
        }

        // 3. Try to extend each existing instance I with the new feature f
        for (const auto& prevInstance : prevInstancesList) {
            // Calculate intersection
            std::vector<const SpatialInstance*> extendedSet = findExtendedSet(
                orderedNRTree, prevInstance, newFeature
            );

            // 4. Create new instances
            for (const auto* neighbor : extendedSet) {
                ColocationInstance newRow = prevInstance;
                newRow.push_back(neighbor);
                newTableRows.push_back(std::move(newRow));
            }
        }

        // Store results
        if (!newTableRows.empty()) {
            result[candidate] = std::move(newTableRows);
        }
        else {
            // --- [DEBUG 4] Xử lý xong nhưng không tạo được instance nào (do không thỏa mãn khoảng cách) ---
            // Đây thường không phải lỗi code, mà do dữ liệu không có co-location ở vị trí này
            std::cout << "[INFO] Candidate ";
            printDebugPattern(candidate);
            std::cout << " processed but NO instances generated (No neighbors satisfy distance).\n";
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
