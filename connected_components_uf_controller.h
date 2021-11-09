#ifndef MST_CONNECTED_COMPONENTS_UF_CONTROLLER_H_
#define MST_CONNECTED_COMPONENTS_UF_CONTROLLER_H_

class ConnectedComponentsUFController {
 public:
  ConnectedComponentsUFController() = delete;
  
  ConnectedComponentsUFController(int num_nodes);
  
  ~ConnectedComponentsUFController();
  
  bool BelongToTheSameConnectedComponent(int u, int v);
  
  void Connect(int u, int v);
  
  bool IsConnected();

 private:
  int GetRepresentative(int u);
  
  const int num_nodes_;
  // representatives_nodes_ids_[i]: representative of the i-th node. Used to track connected components (union-find algorithm)
  int *representatives_nodes_ids_;
  // connected_components_sizes_[i]: the size of the connected component of the i-th node, iff i is its connected component's representative
  int *connected_components_sizes_;
};

#endif  // MST_CONNECTED_COMPONENTS_UF_CONTROLLER_H_
