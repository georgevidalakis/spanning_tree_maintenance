#include "connected_components_dfs_controller.h"

#include <cstdlib>
#include <iostream>
#include <vector>
#include <set>

ConnectedComponentsDFSController::ConnectedComponentsDFSController(const int num_nodes) : num_nodes_(num_nodes), edges_(std::vector<std::set<int>>(num_nodes_)) {}

ConnectedComponentsDFSController::~ConnectedComponentsDFSController() {}

bool ConnectedComponentsDFSController::BelongToTheSameConnectedComponent(const int u, const int v) {
  if (u >= num_nodes_ || v >= num_nodes_) {
    std::cerr << "[ERROR]: ConnectedComponentsDFSController queried to connect non-existent node.\n";
    exit(1);
  }
  bool *visited = (bool*)malloc(num_nodes_ * sizeof(bool));
  for (int node_id = 0; node_id < num_nodes_; ++node_id) {
    visited[node_id] = false;
  }
  DFS(visited, u);
  bool belong_to_the_same_connected_component = visited[v];
  free(visited);
  return belong_to_the_same_connected_component;
}

void ConnectedComponentsDFSController::Connect(const int u, const int v) {
  if (u >= num_nodes_ || v >= num_nodes_) {
    std::cerr << "[ERROR]: ConnectedComponentsDFSController queried to connect non-existent node.\n";
    exit(1);
  }
  edges_[u].insert(v);
  edges_[v].insert(u);
}

bool ConnectedComponentsDFSController::IsConnected() {
  if (num_nodes_ <= 1) {
    return true;
  }
  bool *visited = (bool*)malloc(num_nodes_ * sizeof(bool));
  for (int node_id = 0; node_id < num_nodes_; ++node_id) {
    visited[node_id] = false;
  }
  DFS(visited, 0);
  bool is_connected = true;
  for (int node_id = 1; node_id < num_nodes_; ++node_id) {
    if (!visited[node_id]) {
      is_connected = false;
      break;
    }
  }
  free(visited);
  return is_connected;
}

void ConnectedComponentsDFSController::DFS(bool *visited, int u) {
  visited[u] = true;
  for (int v : edges_[u]) {
    if (!visited[v]) {
      DFS(visited, v);
    }
  }
}
