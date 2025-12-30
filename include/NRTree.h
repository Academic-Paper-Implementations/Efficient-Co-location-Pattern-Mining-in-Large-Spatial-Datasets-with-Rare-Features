#pragma once
#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include <map>
#include <memory>
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
    explicit NRNode(NodeType nodeType) : type(nodeType), instancePtr(nullptr), featureType("") {}

    // Destructor
    ~NRNode() {
        for (auto child : children) delete child;
        children.clear();
    }

    // Rule of 5: Delete copy (manages raw pointers), default move
    NRNode(const NRNode&) = delete;
    NRNode& operator=(const NRNode&) = delete;
    NRNode(NRNode&&) = default;
    NRNode& operator=(NRNode&&) = default;
};

class NRTree {
private:
    std::unique_ptr<NRNode> root;

    // Recursive function to print tree (for debugging purposes)
    void printRecursive(NRNode* node, int level) const;

public:
    NRTree();
    ~NRTree();

    // Rule of 5: Delete copy (unique_ptr member), default move
    NRTree(const NRTree&) = delete;
    NRTree& operator=(const NRTree&) = delete;
    NRTree(NRTree&&) = default;
    NRTree& operator=(NRTree&&) = default;

    // Build tree from NeighborhoodMgr results
    // Features must be sorted by instance count (ascending) according to the paper
    void build(const NeighborhoodMgr& neighMgr, const std::map<FeatureType, int>& featureCounts);

    // Print tree structure for debugging and verification
    void printTree() const;

    // Getter for root if external processing needed
    const NRNode* getRoot() const { return root.get(); }
};