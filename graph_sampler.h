#ifndef MST_GRAPH_SAMPLER_H_
#define MST_GRAPH_SAMPLER_H_

#include <vector>

#include "general.h"
#include "edge_intervals_uniform_sampler.h"

class GraphSampler {
 public:
  GraphSampler() = delete;
  
  GraphSampler(int t, int d, int n, int m);
  
  void AddEdge(int u, int v);
  
  std::vector<Edge> GetEdgesWithRefreshedIntervals();
  
 private:
  const int t_;
  const int d_;
  const int n_;
  const int m_;
  std::vector<Edge> edges_;
  EdgeIntervalsUniformSampler edge_intervals_uniform_sampler_;
};

#endif  // MST_GRAPH_SAMPLER_H_
