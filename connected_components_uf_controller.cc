#include "connected_components_uf_controller.h"

#include <cstdlib>
#include <iostream>

ConnectedComponentsUFController::ConnectedComponentsUFController(const int num_nodes) : num_nodes_(num_nodes) {
  representatives_nodes_ids_ = (int*)malloc(num_nodes_ * sizeof(int));
  connected_components_sizes_ = (int*)malloc(num_nodes_ * sizeof(int));
  if (num_nodes_ && (representatives_nodes_ids_ == NULL || connected_components_sizes_ == NULL)) {
    std::cerr << "[ERROR]: ConnectedComponentsUFController could not allocate the required memory.\n";
    exit(1);
  }
  for (int nodeId = 0; nodeId < num_nodes_; ++nodeId) {
    representatives_nodes_ids_[nodeId] = nodeId;
    connected_components_sizes_[nodeId] = 1;
  }
}

ConnectedComponentsUFController::~ConnectedComponentsUFController() {
  free(representatives_nodes_ids_);
  free(connected_components_sizes_);
}

bool ConnectedComponentsUFController::BelongToTheSameConnectedComponent(const int u, const int v) {
  if (u >= num_nodes_ || v >= num_nodes_) {
    std::cerr << "[ERROR]: ConnectedComponentsUFController queried to connect non-existent node.\n";
    exit(1);
  }
  return (GetRepresentative(u) == GetRepresentative(v));
}

void ConnectedComponentsUFController::Connect(const int u, const int v) {
  if (u >= num_nodes_ || v >= num_nodes_) {
    std::cerr << "[ERROR]: ConnectedComponentsUFController queried to connect non-existent node.\n";
    exit(1);
  }
  const int r_u = GetRepresentative(u), r_v = GetRepresentative(v);
  if (representatives_nodes_ids_[r_u] != r_u || representatives_nodes_ids_[r_v] != r_v) {
    std::cerr << "[ERROR]: ConnectedComponentsUFController DirectConnect queried to connect non-representative node.\n";
    std::exit(1);
  }
  if (r_u != r_v) {
    if (connected_components_sizes_[r_u] >= connected_components_sizes_[r_v]) {
      representatives_nodes_ids_[r_v] = r_u;
      connected_components_sizes_[r_u] += connected_components_sizes_[r_v];
    }
    else {
      representatives_nodes_ids_[r_u] = r_v;
      connected_components_sizes_[r_v] += connected_components_sizes_[r_u];
    }
  }
}

bool ConnectedComponentsUFController::IsConnected() {
  return (!num_nodes_ || connected_components_sizes_[GetRepresentative(0)] == num_nodes_);
}

int ConnectedComponentsUFController::GetRepresentative(int u) {
  while (representatives_nodes_ids_[u] != u) {
    representatives_nodes_ids_[u] = representatives_nodes_ids_[representatives_nodes_ids_[u]];
    u = representatives_nodes_ids_[u];
  }
  return u;
}
