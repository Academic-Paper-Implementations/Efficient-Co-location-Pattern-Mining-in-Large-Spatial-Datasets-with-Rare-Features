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

void NRTree::build(const NeighborhoodMgr& neighMgr) {
    // 0. Reset cây nếu đã có dữ liệu cũ
    if (root) delete root;
    root = new NRNode(ROOT_NODE);

    // Lấy dữ liệu map thô: unordered_map<FeatureType, vector<OrderedNeigh>>
    const auto& rawMap = neighMgr.getOrderedNeighbors();

    // 1. LEVEL 1: FEATURE NODES
    // Vì rawMap là unordered_map, ta cần lấy keys ra và sort để đảm bảo thứ tự cây
    std::vector<FeatureType> sortedFeatures;
    for (const auto& pair : rawMap) {
        sortedFeatures.push_back(pair.first);
    }
    // Sắp xếp Feature theo thứ tự (A, B, C...) hoặc theo quy tắc của bạn
    std::sort(sortedFeatures.begin(), sortedFeatures.end());

    for (const auto& fType : sortedFeatures) {
        // Tạo nút Feature (VD: Node A)
        NRNode* fNode = new NRNode(FEATURE_NODE);
        fNode->featureType = fType;
        root->children.push_back(fNode);

        // Lấy danh sách các Center Instance thuộc Feature này
        const auto& starList = rawMap.at(fType);

        // 2. LEVEL 2: INSTANCE NODES (Center)
        // starList là vector, ta cũng nên đảm bảo nó được sort theo ID instance
        // (Copy ra vector tạm con trỏ để sort nếu cần, hoặc giả định đầu vào đã sort)
        // Ở đây ta duyệt trực tiếp.

        for (const auto& star : starList) {
            NRNode* centerNode = new NRNode(INSTANCE_NODE);
            centerNode->data = star.center; // Lưu con trỏ đến dữ liệu gốc
            fNode->children.push_back(centerNode);

            // 3. LEVEL 3: NEIGHBOR NODES
            // Dữ liệu hàng xóm đang nằm trong: star.neighbors (là unordered_map)
            // Ta cần gom tất cả hàng xóm lại và sắp xếp để đưa vào cây

            std::vector<const SpatialInstance*> allNeighbors;

            for (const auto& mapEntry : star.neighbors) {
                const auto& vec = mapEntry.second;
                allNeighbors.insert(allNeighbors.end(), vec.begin(), vec.end());
            }

            // Sắp xếp hàng xóm: Ưu tiên theo FeatureType, sau đó đến ID
            std::sort(allNeighbors.begin(), allNeighbors.end(),
                [](const SpatialInstance* a, const SpatialInstance* b) {
                    if (a->type != b->type) return a->type < b->type;
                    return a->id < b->id;
                });

            // Tạo node con cho từng hàng xóm
            for (const auto* neighborPtr : allNeighbors) {
                NRNode* neighNode = new NRNode(NEIGHBOR_NODE);
                neighNode->data = neighborPtr;
                centerNode->children.push_back(neighNode);
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
        std::cout << "ROOT (Null)\n";
    }
    else if (node->type == FEATURE_NODE) {
        std::cout << indent << "+ Feature: " << node->featureType << "\n";
    }
    else if (node->type == INSTANCE_NODE) {
        std::cout << indent << "+ Instance: " << node->data->id
            << " (" << node->data->type << ")\n";
    }
    else if (node->type == NEIGHBOR_NODE) {
        std::cout << indent << "- Neighbor: " << node->data->id
            << " [" << node->data->type << "]\n";
        return; // Neighbor là lá, không in tiếp
    }

    // Đệ quy in các con
    for (auto child : node->children) {
        printRecursive(child, level + 1);
    }
}