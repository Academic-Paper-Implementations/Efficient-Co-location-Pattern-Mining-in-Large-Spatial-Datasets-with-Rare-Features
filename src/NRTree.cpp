#include "NRTree.h"
#include "utils.h"

NRTree::NRTree() {
    root = new NRNode(ROOT_NODE);
}

NRTree::~NRTree() {
    if (root) {
        delete root; // Node destructor will automatically delete children recursively
        root = nullptr;
    }
}

void NRTree::build(const NeighborhoodMgr& neighMgr, const std::map<FeatureType, int>& featureCounts, const std::vector<SpatialInstance>& instances) {
    // 0. Reset tree if old data exists
    if (root) delete root;
    root = new NRNode(ROOT_NODE);

    // Get raw map data: unordered_map<FeatureType, vector<OrderedNeigh>>
    const auto& rawMap = neighMgr.getOrderedNeighbors();

    // 1. LEVEL 1: FEATURE NODES
    // According to paper: Features must be sorted by instance count (ascending order)
    // Same logic as in isOrdered() and featureSort()
    std::vector<FeatureType> sortedFeatures;
    for (const auto& pair : rawMap) {
        sortedFeatures.push_back(pair.first);
    }
    //Sort features by feature count (ascending)
    sortedFeatures = featureSort(sortedFeatures, instances);

    for (const auto& fType : sortedFeatures) {
        // Create feature node (e.g., Node A)
        NRNode* fNode = new NRNode(FEATURE_NODE);
        fNode->featureType = fType;
        root->children.push_back(fNode);

        // Get list of center instances for this feature
        const auto& starList = rawMap.at(fType);

        // 2. LEVEL 2: INSTANCE NODES (Center)
        for (const auto& star : starList) {
            NRNode* centerNode = new NRNode(INSTANCE_NODE);
            centerNode->data = star.center; // Store pointer to original data
            fNode->children.push_back(centerNode);

            // 3. LEVEL 3: FEATURE NODES (for neighbor features)
            // Neighbor data is in: star.neighbors (unordered_map<FeatureType, vector<const SpatialInstance*>>)
            // We need to create FEATURE_NODE for each neighbor feature type, sorted by feature count

            // Get list of neighbor feature types and sort by feature count
            std::vector<FeatureType> neighborFeatureTypes;
            for (const auto& mapEntry : star.neighbors) {
                neighborFeatureTypes.push_back(mapEntry.first);
            }
            // Sort neighbor feature types by feature count (ascending)
            neighborFeatureTypes = featureSort(neighborFeatureTypes, instances);

            // Create FEATURE_NODE for each neighbor feature type
            for (const auto& neighborFeatureType : neighborFeatureTypes) {
                NRNode* neighborFeatureNode = new NRNode(FEATURE_NODE);
                neighborFeatureNode->featureType = neighborFeatureType;
                centerNode->children.push_back(neighborFeatureNode);

                // 4. LEVEL 4: INSTANCE_VECTOR_NODE (single node containing vector of neighbor instances)
                // Get list of instances for this feature type and sort by ID (alphabetical)
                const auto& neighborInstances = star.neighbors.at(neighborFeatureType);

                // Create single INSTANCE_VECTOR_NODE to store vector of neighbor instances
                NRNode* instanceVectorNode = new NRNode(INSTANCE_VECTOR_NODE);
                instanceVectorNode->instanceVector = neighborInstances;  // Store entire vector
                neighborFeatureNode->children.push_back(instanceVectorNode);
            }
        }
    }
}

// --- Support functions for display (Debug) ---

void NRTree::printTree() const {
    std::cout << "\n=== ORDERED NR-TREE STRUCTURE ===\n";
    printRecursive(root, 0);
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
        std::cout << indent << "- Instance: " << node->data->id
            << " [" << node->data->type << "]\n";
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