#ifndef MST_GENERAL_H_
#define MST_GENERAL_H_

#include <vector>

struct Interval {
  int start, end;  // the starting and ending time points of the interval, both included
};

struct Edge {
  int u, v;
  std::vector<Interval> intervals;
};

#endif  // MST_GENERAL_H_
