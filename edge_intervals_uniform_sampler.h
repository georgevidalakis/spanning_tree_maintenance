#ifndef MST_EDGE_INTERVALS_UNIFORM_SAMPLER_H_
#define MST_EDGE_INTERVALS_UNIFORM_SAMPLER_H_

#include <stdint.h>

#include <vector>

#include "general.h"

class EdgeIntervalsUniformSampler {
 public:
  EdgeIntervalsUniformSampler() = delete;
  
  EdgeIntervalsUniformSampler(int t, int d);
  
  ~EdgeIntervalsUniformSampler();
  
  // The selection of the time intervals is based on dp_, so that each valid selection has the same probability
  std::vector<Interval> GetRandEdgeIntervals() const;

 private:
  // rem_t: number of remaining time points, rem_t > 0
  int GetRandIntervalLen(int rem_t) const;
  
  const int t_;  // number of time points
  const int d_;  // max possible length of an edge's time-interval in time points
  int64_t *dp_;  // dp_[i]: how many different time-intervals for a single edge exist if t == i
};

#endif  // MST_EDGE_INTERVALS_UNIFORM_SAMPLER_H_
