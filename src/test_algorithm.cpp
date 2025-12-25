#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <algorithm>
#include <cmath>
#include <iomanip>

// Mock FeatureType as string for testing
using FeatureType = std::string;

// ==========================================
// CODE TO TEST (Pasted as requested)
// ==========================================

// Step 2: Sorting features in ascending order of the quantity of instances
std::vector<FeatureType> featureSort(std::vector<FeatureType>& featureSet, const std::map<FeatureType, int>& featureCounts) {
    // Sort logic based on Algorithm 1 Step 2
    // Ascending order of instance counts
    std::sort(featureSet.begin(), featureSet.end(), 
        [&featureCounts](const FeatureType& a, const FeatureType& b) {
            int countA = 0;
            int countB = 0;
            if (featureCounts.find(a) != featureCounts.end()) countA = featureCounts.at(a);
            if (featureCounts.find(b) != featureCounts.end()) countB = featureCounts.at(b);
            
            // Primary sort key: count (ascending)
            if (countA != countB) {
                return countA < countB;
            }
            // Secondary sort key: lexicographical (for stability/determinism)
            return a < b;
        }
    );
    return featureSet;
}

// Step 3: Calculating delta for the spatial dataset
// Formula: delta = (2 / (m*(m-1))) * Sum_{i<j} (num(f_j) / num(f_i))
double calculateDelta(const std::map<FeatureType, int>& featureCounts) {
    if (featureCounts.size() < 2) {
        return 0.0;
    }

    // 1. Extract counts into a vector for indexed access
    std::vector<double> counts;
    for (const auto& pair : featureCounts) {
        counts.push_back(static_cast<double>(pair.second));
    }

    // 2. Sort counts in ascending order (matching Step 2's logic implication)
    // This ensures that for i < j, counts[i] <= counts[j], so the ratio >= 1.
    std::sort(counts.begin(), counts.end());

    double m = static_cast<double>(counts.size());
    double sumRatios = 0.0;

    // 3. Calculate sum of ratios for all pairs i < j
    for (int i = 0; i < m; ++i) {
        for (int j = i + 1; j < m; ++j) {
            // Avoid division by zero if count is 0 (though rare in mining)
            double denominator = counts[i];
            if (denominator == 0) denominator = 1e-6; 
            
            double ratio = counts[j] / denominator;
            sumRatios += ratio;
        }
    }

    // 4. Calculate final Delta
    // Factor = 2 / (m * (m - 1))
    double factor = 2.0 / (m * (m - 1.0));
    
    return factor * sumRatios;
}

// ==========================================
// TEST RUNNER
// ==========================================

void testFeatureSort() {
    std::cout << "\n--- Testing featureSort ---\n";
    
    // Setup data
    std::map<FeatureType, int> counts;
    counts["A"] = 100;
    counts["B"] = 10;   // Rare
    counts["C"] = 50;   // Medium
    
    std::vector<FeatureType> features = {"A", "B", "C"};
    
    std::cout << "Original order: ";
    for(const auto& f : features) std::cout << f << "(" << counts[f] << ") ";
    std::cout << "\n";
    
    // Run sort
    featureSort(features, counts);
    
    std::cout << "Sorted order:   ";
    for(const auto& f : features) std::cout << f << "(" << counts[f] << ") ";
    std::cout << "\n";
    
    // Verify
    bool passed = (features[0] == "B" && features[1] == "C" && features[2] == "A");
    std::cout << "Result: " << (passed ? "PASSED [Correctly sorted Ascending]" : "FAILED") << "\n";
}

void testCalculateDelta() {
    std::cout << "\n--- Testing calculateDelta ---\n";
    
    // Case 1: Simple Pair
    // A=10, B=20
    // Ratio = 20/10 = 2.0
    // Factor = 2 / (2 * 1) = 1.0
    // Expected Delta = 2.0
    std::map<FeatureType, int> counts1;
    counts1["A"] = 10;
    counts1["B"] = 20;
    
    double delta1 = calculateDelta(counts1);
    std::cout << "Case 1 (10, 20): Expected 2.0 | Got " << delta1 << "\n";
    
    // Case 2: Three items
    // A=10, B=20, C=40
    // Pairs (sorted): (10, 20), (10, 40), (20, 40)
    // Ratios: 20/10=2, 40/10=4, 40/20=2
    // Sum Ratios = 2 + 4 + 2 = 8
    // Factor = 2 / (3 * 2) = 1/3
    // Expected Delta = 8 / 3 = 2.666...
    std::map<FeatureType, int> counts2;
    counts2["A"] = 10;
    counts2["B"] = 20;
    counts2["C"] = 40;
    
    double delta2 = calculateDelta(counts2);
    std::cout << "Case 2 (10, 20, 40): Expected 2.666... | Got " << delta2 << "\n";
}

int main() {
    testFeatureSort();
    testCalculateDelta();
    return 0;
}
