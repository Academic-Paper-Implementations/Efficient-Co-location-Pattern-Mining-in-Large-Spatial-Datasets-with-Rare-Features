#include "NRTree.h"

NRTree::NRTree() : root(std::make_unique<NRNode>(ROOT_NODE)) {
}

NRTree::~NRTree() = default;

void NRTree::build(const NeighborhoodMgr& neighMgr, const std::map<FeatureType, int>& featureCounts) {
    // 0. Reset tree if old data exists
    root = std::make_unique<NRNode>(ROOT_NODE);

    // Get raw map data: unordered_map<FeatureType, vector<OrderedNeigh>>
    const auto& rawMap = neighMgr.getOrderedNeighbors();

    // 1. LEVEL 1: FEATURE NODES
    // According to paper: Features must be sorted by instance count (ascending order)
    // Same logic as in isOrdered() and featureSort()
    std::vector<FeatureType> sortedFeatures;
    sortedFeatures.reserve(rawMap.size());
    for (const auto& pair : rawMap) {
        sortedFeatures.push_back(pair.first);
    }
    // Sort features by instance count (ascending), then lexicographic
    std::sort(sortedFeatures.begin(), sortedFeatures.end(),
        [&featureCounts](const FeatureType& a, const FeatureType& b) {
            int countA = featureCounts.count(a) ? featureCounts.at(a) : 0;
            int countB = featureCounts.count(b) ? featureCounts.at(b) : 0;
            if (countA != countB) return countA < countB;
            return a < b; // Lexicographic tie-breaker
        });

    for (const auto& featureType : sortedFeatures) {
        // Create feature node (e.g., Node A)
        NRNode* featureNode = new NRNode(FEATURE_NODE);
        featureNode->featureType = featureType;
        root->children.push_back(featureNode);

        // Get list of center instances for this feature
        const auto& starList = rawMap.at(featureType);

        // 2. LEVEL 2: INSTANCE NODES (Center)
        // Sort instance centers by ID for consistent ordering
        std::vector<OrderedNeigh> sortedStarList = starList;
        std::sort(sortedStarList.begin(), sortedStarList.end(),
            [](const OrderedNeigh& a, const OrderedNeigh& b) {
                return a.center->id < b.center->id;
            });
		/////////////////////////////////////////////////////////
        // TODO: edit sort by other way, this way is wrong.
		/////////////////////////////////////////////////////////

        for (const auto& star : sortedStarList) {
            NRNode* centerNode = new NRNode(INSTANCE_NODE);
            centerNode->instancePtr = star.center; // Store pointer to original data
            featureNode->children.push_back(centerNode);

            // 3. LEVEL 3: FEATURE NODES (for neighbor features)
            // Neighbor data is in: star.neighbors (unordered_map<FeatureType, vector<const SpatialInstance*>>)
            // We need to create FEATURE_NODE for each neighbor feature type, sorted by feature count

            // Get list of neighbor feature types and sort by feature count
            std::vector<FeatureType> neighborFeatureTypes;
            neighborFeatureTypes.reserve(star.neighbors.size());
            for (const auto& mapEntry : star.neighbors) {
                neighborFeatureTypes.push_back(mapEntry.first);
            }

            // Sort neighbor feature types by feature count (ascending), then lexicographic
            std::sort(neighborFeatureTypes.begin(), neighborFeatureTypes.end(),
                [&featureCounts](const FeatureType& a, const FeatureType& b) {
                    int countA = featureCounts.count(a) ? featureCounts.at(a) : 0;
                    int countB = featureCounts.count(b) ? featureCounts.at(b) : 0;
                    if (countA != countB) return countA < countB;
                    return a < b; // Lexicographic tie-breaker
                });
            /////////////////////////////////////////////////////////
            // TODO: change this sort function to utils.h.
            /////////////////////////////////////////////////////////

            // Create FEATURE_NODE for each neighbor feature type
            for (const auto& neighborFeatureType : neighborFeatureTypes) {
                NRNode* neighborFeatureNode = new NRNode(FEATURE_NODE);
                neighborFeatureNode->featureType = neighborFeatureType;
                centerNode->children.push_back(neighborFeatureNode);

                // 4. LEVEL 4: INSTANCE_VECTOR_NODE (single node containing vector of neighbor instances)
                // Get list of instances for this feature type and sort by ID (alphabetical)
                const auto& neighborInstances = star.neighbors.at(neighborFeatureType);
                std::vector<const SpatialInstance*> sortedNeighborInstances = neighborInstances;

                // Sort by ID (alphabetical order)
                std::sort(sortedNeighborInstances.begin(), sortedNeighborInstances.end(),
                [](const SpatialInstance* a, const SpatialInstance* b) {
                    if (a->type != b->type) return a->type < b->type;
                    return a->id < b->id;
                });
                /////////////////////////////////////////////////////////
                // TODO: edit sort by other way, this way is wrong.
                /////////////////////////////////////////////////////////

                // Create single INSTANCE_VECTOR_NODE to store vector of neighbor instances
                NRNode* instanceVectorNode = new NRNode(INSTANCE_VECTOR_NODE);
                instanceVectorNode->instanceVector = sortedNeighborInstances;  // Store entire vector
                neighborFeatureNode->children.push_back(instanceVectorNode);
            }
        }
    }
}

// --- Support functions for display (Debug) ---

void NRTree::printTree() const {
    std::cout << "\n=== ORDERED NR-TREE STRUCTURE ===\n";
    printRecursive(root.get(), 0);
    std::cout << "=================================\n";
}

void NRTree::printRecursive(NRNode* node, int level) const {
    // Indent according to level
    std::string indent = "";
    for (int i = 0; i < level; i++) indent += "  | ";

    if (node->type == ROOT_NODE) {
        std::cout << "ROOT\n";
    }
    else if (node->type == FEATURE_NODE) {
        std::cout << indent << "+ Feature: " << node->featureType << "\n";
    }
    else if (node->type == INSTANCE_NODE) {
        std::cout << indent << "- Instance: " << node->instancePtr->id
            << " [" << node->instancePtr->type << "]\n";
    }
    else if (node->type == INSTANCE_VECTOR_NODE) {
        std::cout << indent << "- Instance Vector (" << node->instanceVector.size() << " instances): [";
        bool first = true;
        for (const auto* inst : node->instanceVector) {
            if (!first) std::cout << ", ";
            std::cout << inst->id << "[" << inst->type << "]";
            first = false;
    }
        std::cout << "]\n";
        return; // INSTANCE_VECTOR_NODE is a leaf, no children
    }

    // New tree structure:
    // Level 1: FEATURE_NODE (center feature, sorted by feature count)
    // Level 2: INSTANCE_NODE (center instance, sorted by ID)
    // Level 3: FEATURE_NODE (neighbor feature, sorted by feature count)
    // Level 4: INSTANCE_VECTOR_NODE (vector of neighbor instances, sorted alphabetically by ID)
    // INSTANCE_VECTOR_NODE is a leaf, no children

    // Recursively print children
    for (auto child : node->children) {
        printRecursive(child, level + 1);
    }
}