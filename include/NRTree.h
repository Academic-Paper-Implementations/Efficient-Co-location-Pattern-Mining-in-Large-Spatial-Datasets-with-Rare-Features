#pragma once
#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include <map>
#include "neighborhood_mgr.h" // Để dùng struct OrderedNeigh và FeatureType
#include "types.h"

// --- [QUAN TRỌNG] THÊM DÒNG NÀY ---
class NeighborhoodMgr;
// ---------------------------------

// Loại node trong cây để dễ quản lý
enum NodeType {
    ROOT_NODE,
    FEATURE_NODE,
    INSTANCE_NODE,
    INSTANCE_VECTOR_NODE,  // Node để lưu vector các instances (dùng ở Level 4)
    NEIGHBOR_NODE
};

// Cấu trúc một nút trong cây NR-Tree
struct NRNode {
    NodeType type;

    // Dữ liệu tùy thuộc vào loại node
    FeatureType featureType;        // Dùng nếu là FEATURE_NODE
    const SpatialInstance* data;    // Dùng nếu là INSTANCE_NODE hoặc NEIGHBOR_NODE
    std::vector<const SpatialInstance*> instanceVector;  // Dùng nếu là INSTANCE_VECTOR_NODE

    // Danh sách con
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

    // Hàm đệ quy để in cây (cho mục đích debug)
    void printRecursive(NRNode* node, int level) const;

public:
    inline NRTree();
    ~NRTree();

    // Hàm quan trọng nhất: Xây dựng cây từ kết quả của NeighborhoodMgr
    // Theo paper: features phải được sắp xếp theo số lượng instance (ascending)
    void build(const NeighborhoodMgr& neighMgr, const std::map<FeatureType, int>& featureCounts);

    // Hàm in cây ra màn hình để kiểm tra
    void printTree() const;

    // Getter root nếu cần xử lý bên ngoài
    const NRNode* getRoot() const { return root; }
};