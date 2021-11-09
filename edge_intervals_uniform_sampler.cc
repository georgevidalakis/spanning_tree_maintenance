#include "edge_intervals_uniform_sampler.h"

#include <stdint.h>

#include <cstdlib>
#include <iostream>
#include <vector>

#include "general.h"

EdgeIntervalsUniformSampler::EdgeIntervalsUniformSampler(const int t, const int d) : t_(t), d_(d) {
  dp_ = (int64_t*)malloc((t_ + 1) * sizeof(int64_t));
  if (dp_ == NULL) {
    std::cerr << "[ERROR]: EdgeIntervalsUniformSampler could not allocate the required memory.\n";
    exit(1);
  }
  // If t_ == 0, there is only one way to select an edge's time-intervals, do not select any
  dp_[0] = 1;
  for (int cur_t = 1; cur_t <= t_; ++cur_t) {  // Solve the subproblem with cur_t time points
    int64_t cur_dp = 0;
    // If we do not cover the first time point
    cur_dp += dp_[cur_t - 1];
    const int max_l = std::min(d_, cur_t);
    for (int l = 1; l <= max_l; ++l) {  // If we cover the first l time points with an interval
      cur_dp += dp_[cur_t - l];
    }
    dp_[cur_t] = cur_dp;
  }
}

EdgeIntervalsUniformSampler::~EdgeIntervalsUniformSampler() {
  free(dp_);
}

// The selection of the time intervals is based on dp_, so that each valid selection has the same probability
std::vector<Interval> EdgeIntervalsUniformSampler::GetRandEdgeIntervals() const {
  // cur_t indicates the current time point. It is updated after every choice
  int cur_t = 0;
  std::vector<Interval> edge_intervals;
  while (cur_t < t_) {
    const int l = GetRandIntervalLen(t_ - cur_t);  // l is the length of the selected time interval starting at time point cur_t. If l == 0, then time point cur_t is not covered
    if (!l) {
      ++cur_t;
    }
    else {
      Interval cur_interval;
      cur_interval.start = cur_t;
      cur_interval.end = cur_t + l - 1;
      edge_intervals.push_back(cur_interval);
      cur_t += l;
    }
  }
  return edge_intervals;
}

// rem_t: number of remaining time points, rem_t > 0
int EdgeIntervalsUniformSampler::GetRandIntervalLen(int rem_t) const {
  const int max_l = std::min(rem_t, d_);
  int64_t sum = 0;
  sum += dp_[rem_t - 1];  // If l == 0, we just leave one time point uncovered
  for (int l = 1; l <= max_l; ++l) {
    sum += dp_[rem_t - l];  // We cover l time points
  }
  double *probs = (double*)malloc((max_l + 1) * sizeof(double));  // probs[i]: the probability that l == i is selected
  if (probs == NULL) {
    std::cerr << "[ERROR]: EdgeIntervalsUniformSampler's GetRandIntervalLen could not allocate the required memory.\n";
    exit(1);
  }
  // probabilities are proportional to the ways to select the time-intervals after l's selection (in order for every valid edge's interval selection to have equal probability)
  probs[0] = (double)dp_[rem_t - 1] / sum;
  for (int l = 1; l <= max_l; ++l) {
    probs[l] = (double)dp_[rem_t - l] / sum;
  }
  double rand_double = (double)rand() / RAND_MAX;
  int selected_l = max_l;
  for (int l = 0; l < max_l; ++l) {
    rand_double -= probs[l];
    if (rand_double <= 0.0) {
      selected_l = l;
      break;
    }
  }
  free(probs);
  return selected_l;
}
