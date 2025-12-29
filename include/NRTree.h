#pragma once
#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include <map>
#include "neighborhood_mgr.h" // To use struct OrderedNeigh and FeatureType
#include "types.h"

// --- [IMPORTANT] FORWARD DECLARATION ---
class NeighborhoodMgr;
// ---------------------------------

// Node types in the NR-Tree for structured organization
enum NodeType {
    ROOT_NODE,              // Root of the tree
    FEATURE_NODE,           // Node representing a feature type
    INSTANCE_NODE,          // Node representing a spatial instance
    INSTANCE_VECTOR_NODE    // Leaf node containing a vector of instances
};

// Structure of a node in the NR-Tree
struct NRNode {
    NodeType type;

    // Data depends on node type
    FeatureType featureType;                            // Used if FEATURE_NODE
    const SpatialInstance* instancePtr;                 // Used if INSTANCE_NODE
    std::vector<const SpatialInstance*> instanceVector; // Used if INSTANCE_VECTOR_NODE

    // List of children
    std::vector<NRNode*> children;

    // Constructor
    NRNode(NodeType nodeType) : type(nodeType), instancePtr(nullptr), featureType("") {}

    ~NRNode() {
        for (auto child : children) delete child;
        children.clear();
    }
};

class NRTree {
private:
    NRNode* root;

    // Recursive function to print tree (for debugging purposes)
    void printRecursive(NRNode* node, int level) const;

public:
    NRTree();
    ~NRTree();

    // Build tree from NeighborhoodMgr results
    // Features must be sorted by instance count (ascending) according to the paper
    void build(const NeighborhoodMgr& neighMgr, const std::map<FeatureType, int>& featureCounts);

    // Print tree structure for debugging and verification
    void printTree() const;

    // Getter for root if external processing needed
    const NRNode* getRoot() const { return root; }
};