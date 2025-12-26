#include "NRTree.h"

NRTree::NRTree() {
    root = new NRNode(ROOT_NODE);
}

NRTree::~NRTree() {
    if (root) {
        delete root; // Destructor của Node sẽ tự động delete các con đệ quy
        root = nullptr;
    }
}

void NRTree::build(const NeighborhoodMgr& neighMgr, const std::map<FeatureType, int>& featureCounts) {
    // 0. Reset cây nếu đã có dữ liệu cũ
    if (root) delete root;
    root = new NRNode(ROOT_NODE);

    // Lấy dữ liệu map thô: unordered_map<FeatureType, vector<OrderedNeigh>>
    const auto& rawMap = neighMgr.getOrderedNeighbors();

    // 1. LEVEL 1: FEATURE NODES
    // Theo paper: Features phải được sắp xếp theo số lượng instance (ascending order)
    // Giống với logic trong isOrdered() và featureSort()
    std::vector<FeatureType> sortedFeatures;
    for (const auto& pair : rawMap) {
        sortedFeatures.push_back(pair.first);
    }
    // Sắp xếp Feature theo số lượng instance (ascending), sau đó theo lexicographic
    std::sort(sortedFeatures.begin(), sortedFeatures.end(),
        [&featureCounts](const FeatureType& a, const FeatureType& b) {
            int countA = featureCounts.count(a) ? featureCounts.at(a) : 0;
            int countB = featureCounts.count(b) ? featureCounts.at(b) : 0;
            if (countA != countB) return countA < countB;
            return a < b; // Lexicographic tie-breaker
        });

    for (const auto& fType : sortedFeatures) {
        // Tạo nút Feature (VD: Node A)
        NRNode* fNode = new NRNode(FEATURE_NODE);
        fNode->featureType = fType;
        root->children.push_back(fNode);

        // Lấy danh sách các Center Instance thuộc Feature này
        const auto& starList = rawMap.at(fType);

        // 2. LEVEL 2: INSTANCE NODES (Center)
        // Sắp xếp các instance centers theo ID để đảm bảo thứ tự nhất quán
        std::vector<OrderedNeigh> sortedStarList = starList;
        std::sort(sortedStarList.begin(), sortedStarList.end(),
            [](const OrderedNeigh& a, const OrderedNeigh& b) {
                return a.center->id < b.center->id;
            });

        for (const auto& star : sortedStarList) {
            NRNode* centerNode = new NRNode(INSTANCE_NODE);
            centerNode->data = star.center; // Lưu con trỏ đến dữ liệu gốc
            fNode->children.push_back(centerNode);

            // 3. LEVEL 3: FEATURE NODES (cho neighbor features)
            // Dữ liệu hàng xóm đang nằm trong: star.neighbors (là unordered_map<FeatureType, vector<const SpatialInstance*>>)
            // Ta cần tạo FEATURE_NODE cho mỗi feature type của neighbor, sắp xếp theo feature count

            // Lấy danh sách các feature types của neighbors và sắp xếp theo feature count
            std::vector<FeatureType> neighborFeatureTypes;
            for (const auto& mapEntry : star.neighbors) {
                neighborFeatureTypes.push_back(mapEntry.first);
            }

            // Sắp xếp neighbor feature types theo feature count (ascending), sau đó lexicographic
            std::sort(neighborFeatureTypes.begin(), neighborFeatureTypes.end(),
                [&featureCounts](const FeatureType& a, const FeatureType& b) {
                    int countA = featureCounts.count(a) ? featureCounts.at(a) : 0;
                    int countB = featureCounts.count(b) ? featureCounts.at(b) : 0;
                    if (countA != countB) return countA < countB;
                    return a < b; // Lexicographic tie-breaker
                });

            // Tạo FEATURE_NODE cho mỗi neighbor feature type
            for (const auto& neighborFeatureType : neighborFeatureTypes) {
                NRNode* neighborFeatureNode = new NRNode(FEATURE_NODE);
                neighborFeatureNode->featureType = neighborFeatureType;
                centerNode->children.push_back(neighborFeatureNode);

                // 4. LEVEL 4: INSTANCE_VECTOR_NODE (một node duy nhất chứa vector các neighbor instances)
                // Lấy danh sách instances của feature type này và sắp xếp theo ID (alphabetical)
                const auto& neighborInstances = star.neighbors.at(neighborFeatureType);
                std::vector<const SpatialInstance*> sortedNeighborInstances = neighborInstances;

                // Sắp xếp theo ID (alphabetical order)
                std::sort(sortedNeighborInstances.begin(), sortedNeighborInstances.end(),
                [](const SpatialInstance* a, const SpatialInstance* b) {
                    if (a->type != b->type) return a->type < b->type;
                    return a->id < b->id;
                });

                // Tạo một INSTANCE_VECTOR_NODE duy nhất để lưu vector các neighbor instances
                NRNode* instanceVectorNode = new NRNode(INSTANCE_VECTOR_NODE);
                instanceVectorNode->instanceVector = sortedNeighborInstances;  // Lưu toàn bộ vector
                neighborFeatureNode->children.push_back(instanceVectorNode);
            }
        }
    }
}

// --- Phần hỗ trợ hiển thị (Debug) ---

void NRTree::printTree() const {
    std::cout << "\n=== ORDERED NR-TREE STRUCTURE ===\n";
    printRecursive(root, 0);
    std::cout << "=================================\n";
}

void NRTree::printRecursive(NRNode* node, int level) const {
    // Thụt đầu dòng theo cấp độ
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
        return; // INSTANCE_VECTOR_NODE là lá, không có con
    }

    // Cấu trúc cây mới:
    // Level 1: FEATURE_NODE (center feature, sorted by feature count)
    // Level 2: INSTANCE_NODE (center instance, sorted by ID)
    // Level 3: FEATURE_NODE (neighbor feature, sorted by feature count)
    // Level 4: INSTANCE_VECTOR_NODE (vector các neighbor instances, sorted alphabetically by ID)
    // INSTANCE_VECTOR_NODE là lá, không có con

    // Đệ quy in các con
    for (auto child : node->children) {
        printRecursive(child, level + 1);
    }
}