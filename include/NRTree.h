#include <vector>
#include <algorithm>


// Create class N-ary tree
class NRTree {
private:
    // CẤU TRÚC LÕI (Core Data Structure)
    // Index của vector cha chính là ID của Instance cha (ví dụ: ID của D.1)
    // Vector con chứa danh sách ID của các Instance hàng xóm (ví dụ: A.1, B.3, C.1)
    std::vector<std::vector<int>> adj_list;

public:
    // Khởi tạo kích thước cây dựa trên tổng số lượng instance tối đa của dữ liệu
    NRTree(int max_instance_id) {
        // +1 để có thể truy cập bằng index ID trực tiếp
        adj_list.resize(max_instance_id + 1);
    }

    // Hàm thêm mối quan hệ cha -> con (Tương ứng các đường nối trong ảnh)
    // parent_id: Ví dụ ID của D.1
    // child_id: Ví dụ ID của A.1
    void insert_relation(int parent_id, int child_id) {
        if (parent_id >= adj_list.size()) {
            // Resize nếu ID vượt quá dự tính (An toàn)
            adj_list.resize(parent_id + 100);
        }
        adj_list[parent_id].push_back(child_id);
    }

    // Hàm lấy danh sách con (Dùng cho bước sinh bảng Instance Table)
    // Tương đương với việc nhìn vào D.1 và thấy A.1, B.3, C.1
    const std::vector<int>& get_children(int parent_id) const {
        if (parent_id < adj_list.size()) {
            return adj_list[parent_id];
        }
        static const std::vector<int> empty;
        return empty;
    }

    // Hàm tiện ích: Sắp xếp lại các con (Đảm bảo tính "Ordered" của cây)
    // Dùng sau khi add hết neighbors
    void sort_children() {
        for (auto& children : adj_list) {
            // Trong thực tế, bạn cần sort theo Feature Type, sau đó đến Instance ID
            // Ở đây sort theo ID để đơn giản hóa demo
            std::sort(children.begin(), children.end());
        }
    }
};
