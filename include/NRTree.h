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

// Node type in tree for easy management
enum NodeType {
    ROOT_NODE,
    FEATURE_NODE,
    INSTANCE_NODE,
    INSTANCE_VECTOR_NODE
};

// Structure of a node in the NR-Tree
struct NRNode {
    NodeType type;

    // Data depends on node type
    FeatureType featureType;        // Used if FEATURE_NODE
    const SpatialInstance* data;    // Used if INSTANCE_NODE or NEIGHBOR_NODE
    std::vector<const SpatialInstance*> instanceVector;  // Used if INSTANCE_VECTOR_NODE

    // List of children
    std::vector<NRNode*> children;

    // Constructor helper
    NRNode(NodeType t) : type(t), data(nullptr), featureType("") {}

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

    // Most important function: Build tree from NeighborhoodMgr results
    // According to paper: features must be sorted by instance count (ascending)
    void build(const NeighborhoodMgr& neighMgr, const std::map<FeatureType, int>& featureCounts);

    // Function to print tree to screen for verification
    void printTree() const;

    // Getter for root if external processing needed
    const NRNode* getRoot() const { return root; }
};