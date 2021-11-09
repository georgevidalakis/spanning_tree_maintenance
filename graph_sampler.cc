#include "graph_sampler.h"

#include <cstdlib>
#include <iostream>
#include <vector>

#include "general.h"
#include "edge_intervals_uniform_sampler.h"

GraphSampler::GraphSampler(const int t, const int d, const int n, const int m) : t_(t), d_(d), n_(n), m_(m), edge_intervals_uniform_sampler_(EdgeIntervalsUniformSampler(t_, d_)) {}

void GraphSampler::AddEdge(const int u, const int v) {
  if ((int)edges_.size() == m_) {
    std::cerr << "[BUG]: GraphSampler tried to surpass m edges.\n";
    exit(1);
  }
  Edge edge;
  edge.u = u;
  edge.v = v;
  edges_.push_back(edge);
}

std::vector<Edge> GraphSampler::GetEdgesWithRefreshedIntervals() {
  if ((int)edges_.size() < m_) {
    std::cerr << "[BUG]: GraphSampler tried to get edges before they bacame m.\n";
    exit(1);
  }
  for (int edge_id = 0; edge_id < m_; ++edge_id) {
    edges_[edge_id].intervals = edge_intervals_uniform_sampler_.GetRandEdgeIntervals();
  }
  return edges_;
}
