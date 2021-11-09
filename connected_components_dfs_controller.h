#ifndef MST_CONNECTED_COMPONENTS_DFS_CONTROLLER_H_
#define MST_CONNECTED_COMPONENTS_DFS_CONTROLLER_H_

#include <vector>
#include <set>

class ConnectedComponentsDFSController {
 public:
  ConnectedComponentsDFSController() = delete;
  
  ConnectedComponentsDFSController(int num_nodes);
  
  ~ConnectedComponentsDFSController();
  
  bool BelongToTheSameConnectedComponent(int u, int v);
  
  void Connect(int u, int v);
  
  bool IsConnected();

 private:
  void DFS(bool *visited, int u);
  
  const int num_nodes_;
  // edges_[i]: the nodes with whom the i-th node is connected
  std::vector<std::set<int>> edges_;
};

#endif  // MST_CONNECTED_COMPONENTS_DFS_CONTROLLER_H_
