#include "solvers.h"

#include <vector>
#include <set>
#include <cstddef>
#include <ctime>
#include <cstdlib>
#include <utility>
#include <algorithm>

#include "connected_components_uf_controller.h"
#include "connected_components_dfs_controller.h"
#include "subset_iterator.h"

void InitRand() {
  srand(time(NULL));
  // In some systems, for small changes of time(NULL), the first few values of rand() are close to each other
  // To increase randomness we ignore them
  for (int i = 0; i < 100; ++i) {
    rand();
  }
}

namespace {

// edge_pair_i: a (interval end, edge id) pair
// This function is used by the greedy algorithm to sort the edges so that those whose current interval ends later are considered first
// The second part of each pair (edge id) is ignored, in order to allow randomness (e.g. if the vector had been shuffled before sorting)
bool CompareEdgesForGreedy(std::pair<int, int> edge_pair_0, std::pair<int, int> edge_pair_1) {
  return edge_pair_0.first > edge_pair_1.first;
}

std::set<int> GetRelevantEdgesIdsSubset(
    std::vector<Edge> &edges, std::set<int> &prev_time_point_edges_ids, std::vector<std::vector<Interval>::iterator> &cur_intervals_its, const int cur_t) {
  std::set<int> relevant_edges_ids;
  for (int edge_id : prev_time_point_edges_ids) {
    if (cur_intervals_its[edge_id] != edges[edge_id].intervals.end() && cur_intervals_its[edge_id]->start <= cur_t - 1) {
      relevant_edges_ids.insert(edge_id);
    }
  }
  return relevant_edges_ids;
}

int CountNewEdges(std::set<int> &relevant_edges_ids, std::set<int> &covering_combination) {
  int num_new_edges = 0;
  for (int edge_id : covering_combination) {
    if (relevant_edges_ids.find(edge_id) == relevant_edges_ids.end()) {
      ++num_new_edges;
    }
  }
  return num_new_edges;
}

bool CouldBeGeneralStepCheck(std::set<int> &relevant_edges_ids, std::set<int> &covering_combination) {
  for (int edge_id : relevant_edges_ids) {
    if (covering_combination.find(edge_id) == covering_combination.end()) {
      return false;
    }
  }
  return true;
}

}  // namespace

// t: number of time points
// n: number of nodes
// m: number of edges
// edges: a vector with the graph's edges
// Check if an instance is feasible
bool IsFeasibleCheck(const int t, const int n, const int m, std::vector<Edge> &edges) {
  std::vector<ConnectedComponentsDFSController*> cc_controllers_ptrs;
  cc_controllers_ptrs.reserve(t);
  for (int cur_t = 0; cur_t < t; ++cur_t) {
    cc_controllers_ptrs.push_back(new ConnectedComponentsDFSController(n));
  }
  for (int edge_id = 0; edge_id < m; ++edge_id) {
    for (Interval interval : edges[edge_id].intervals) {
      for (int cur_t = interval.start; cur_t <= interval.end; ++cur_t) {
        cc_controllers_ptrs[cur_t]->Connect(edges[edge_id].u, edges[edge_id].v);
      }
    }
  }
  bool is_feasible = true;
  for (ConnectedComponentsDFSController* cc_controller_ptr : cc_controllers_ptrs) {
    if (!cc_controller_ptr->IsConnected()) {
      is_feasible = false;
      break;
    }
  }
  for (ConnectedComponentsDFSController* cc_controller_ptr : cc_controllers_ptrs) {
    delete cc_controller_ptr;
  }
  return is_feasible;
}

// t: number of time points
// n: number of nodes
// m: number of edges
// edges: a vector with the graph's edges
// Get the cost (number of selected time-intervals) of the greedy algorithm's solution, or -1 if the instance is infeasible
int SolveGreedy(const int t, const int n, const int m, std::vector<Edge> &edges) {
  int greedy_cost = 0;  // The number of different edge time-intervals selected by the greedy algorithm
  // cur_intervals_its[i]: an iterator to the i-th edge's first time-interval that ends at or after time point cur_t.
  // If there is no such time-interval then cur_intervals_its[i] == edges[i].intervals.end()
  std::vector<std::vector<Interval>::iterator> cur_intervals_its;
  for (int edge_id = 0; edge_id < m; ++edge_id) {
    cur_intervals_its.push_back(edges[edge_id].intervals.begin());
  }
  std::vector<int> selected_edges_ids;  // The ids of the currently selected edges. It is updated every time point
  for (int cur_t = 0; cur_t < t; ++cur_t) {  // cur_t is the current time point
    // Keep already selected edges whose previous time-interval includes the current time point
    for (unsigned int i = 0; i < selected_edges_ids.size(); ++i) {
      const int selected_edge_id = selected_edges_ids[i];
      if (cur_intervals_its[selected_edge_id] == edges[selected_edge_id].intervals.end() || cur_intervals_its[selected_edge_id]->end < cur_t) {
        selected_edges_ids[i] = selected_edges_ids.back();
        selected_edges_ids.pop_back();
        --i;
      }
    }
    // Update cur_intervals_its
    for (int edge_id = 0; edge_id < m; ++edge_id) {
      if (cur_intervals_its[edge_id] != edges[edge_id].intervals.end() && cur_intervals_its[edge_id]->end < cur_t) {
        ++cur_intervals_its[edge_id];
      }
    }
    // Initialize the connected components
    ConnectedComponentsUFController cc_controller(n);
    for (int selected_edge_id : selected_edges_ids) {
      cc_controller.Connect(edges[selected_edge_id].u, edges[selected_edge_id].v);
    }
    // Find which non selected edges are available at the current time point and sort them by the ending time of their current interval in descending order
    std::vector<std::pair<int, int>> candidate_edges;  // (interval end, edge id) pairs
    for (int edge_id = 0; edge_id < m; ++edge_id) {
      if (cur_intervals_its[edge_id] != edges[edge_id].intervals.end() && cur_intervals_its[edge_id]->start <= cur_t) {
        candidate_edges.push_back(std::make_pair(cur_intervals_its[edge_id]->end, edge_id));
      }
    }
    random_shuffle(candidate_edges.begin(), candidate_edges.end());
    sort(candidate_edges.begin(), candidate_edges.end(), CompareEdgesForGreedy);
    // Consider the candidate edges one by one and select those whose nodes belong to different connected components
    // When an edge is selected the connected components of its nodes are merged
    for (std::pair<int, int> candidate_edge : candidate_edges) {
      const int edge_id = candidate_edge.second;
      if (!cc_controller.BelongToTheSameConnectedComponent(edges[edge_id].u, edges[edge_id].v)) {
        selected_edges_ids.push_back(edge_id);
        cc_controller.Connect(edges[edge_id].u, edges[edge_id].v);
        // A new time-interval has been selected by the greedy algorithm
        ++greedy_cost;
      }
    }
    // Check whether the instance is feasible for the current time point
    if (!cc_controller.IsConnected()) {
      // cout << "The instance is infeasible for the " << cur_t << "-th time point (according to the greedy algorithm).\n";
      return -1;
    }
  }
  return greedy_cost;
}

// t: number of time points
// n: number of nodes
// m: number of edges
// edges: a vector with the graph's edges
// Get the cost (number of selected time-intervals) of an optimal solution, or -1 if the instance is infeasible
int SolveOptimal(const int t, const int n, const int m, std::vector<Edge> &edges) {
  // Try every possible interval selection and keep the lowest cost (if there is any feasible solution)
  std::vector<std::pair<int, std::vector<Interval>::iterator>> candidate_intervals;  // (edge id, edge interval it) pairs
  for (int edge_id = 0; edge_id < m; ++edge_id) {
    for (std::vector<Interval>::iterator edge_interval_it = edges[edge_id].intervals.begin(); edge_interval_it != edges[edge_id].intervals.end(); ++edge_interval_it) {
      candidate_intervals.push_back(std::make_pair(edge_id, edge_interval_it));
    }
  }
  int optimal_cost = -1;
  SubsetIterator<std::pair<int, std::vector<Interval>::iterator>> selectedAvailableEdgesSubsetIterator(candidate_intervals);
  while (!selectedAvailableEdgesSubsetIterator.HasReachedEndOfIterator()) {
    std::set<std::pair<int, std::vector<Interval>::iterator>> selected_intervals = selectedAvailableEdgesSubsetIterator.GetSubset();  // (edge id, edge interval it) pairs
    // Check if the selected intervals consist a feasible solution
    bool is_feasible = true;
    for (int cur_t = 0; cur_t < t; ++cur_t) {
      ConnectedComponentsDFSController cc_controller(n);
      for (std::pair<int, std::vector<Interval>::iterator> selected_interval : selected_intervals) {
        const int edge_id = selected_interval.first;
        const Interval interval = *(selected_interval.second);
        if (interval.start <= cur_t && cur_t <= interval.end) {
          cc_controller.Connect(edges[edge_id].u, edges[edge_id].v);
        }
      }
      if (!cc_controller.IsConnected()) {
        is_feasible = false;
        break;
      }
    }
    if (is_feasible && (optimal_cost == -1 || (int)selected_intervals.size() < optimal_cost)) {
      optimal_cost = selected_intervals.size();
    }
  }
  return optimal_cost;
}

bool CouldBeOptimalStepCheck(
    const int n,
    std::vector<Edge> &edges,
    std::set<int> &relevant_edges_ids,
    std::set<int> &covering_combination,
    std::vector<int> &available_edges_ids,
    std::vector<std::vector<Interval>::iterator> &cur_intervals_its) {
  return true;
}

bool CouldBeGreedyStepCheck(
    const int n,
    std::vector<Edge> &edges,
    std::set<int> &relevant_edges_ids,
    std::set<int> &covering_combination,
    std::vector<int> &available_edges_ids,
    std::vector<std::vector<Interval>::iterator> &cur_intervals_its) {
  // Create connected components before the selection of new edges
  ConnectedComponentsUFController cc_controller(n);
  for (int edge_id : relevant_edges_ids) {
    cc_controller.Connect(edges[edge_id].u, edges[edge_id].v);
  }
  // Order the available edges by their current interval end time-point in descending order, but discard those already used
  std::vector<std::pair<int, int>> available_edges;  // (interval end, edge id) pairs
  for (int edge_id : available_edges_ids) {
    if (relevant_edges_ids.find(edge_id) == relevant_edges_ids.end()) {
      available_edges.push_back(std::make_pair(cur_intervals_its[edge_id]->end, edge_id));
    }
  }
  sort(available_edges.begin(), available_edges.end(), CompareEdgesForGreedy);
  // Order the selected edges by their current interval end time-point in descending order, but discard those already used
  std::vector<std::pair<int, int>> selected_edges;  // (interval end, edge id) pairs
  for (int edge_id : covering_combination) {
    if (relevant_edges_ids.find(edge_id) == relevant_edges_ids.end()) {
      selected_edges.push_back(std::make_pair(cur_intervals_its[edge_id]->end, edge_id));
    }
  }
  sort(selected_edges.begin(), selected_edges.end(), CompareEdgesForGreedy);
  // The step could be greedy if:
  // a) The selected edges don't cause any cicle.
  // b) The selected edges used as connections in order (after the sorting above) don't have earlier ending time-points than any not causing cicles available edge
  int available_edge_idx = 0;
  for (std::pair<int, int> selected_edge : selected_edges) {
    int selected_edge_id = selected_edge.second;
    if (cc_controller.BelongToTheSameConnectedComponent(edges[selected_edge_id].u, edges[selected_edge_id].v)) {
      return false;
    }
    while (true) {
      int available_edge_id = available_edges[available_edge_idx].second;
      if (cc_controller.BelongToTheSameConnectedComponent(edges[available_edge_id].u, edges[available_edge_id].v)) {
        ++available_edge_idx;
        continue;
      }
      if (selected_edge.first < available_edges[available_edge_idx].first) {
        return false;
      }
      break;
    }
    cc_controller.Connect(edges[selected_edge_id].u, edges[selected_edge_id].v);
  }
  return true;
}

// t: number of time points
// n: number of nodes
// m: number of edges
// edges: a vector with the graph's edges
// could_be_specific_step_check: function that checks if a combination of edges selected for two consequitive time-points is possible for an algorithm (e.g. optimal, greedy)
// compute_best: Specifies whether to compute the best or the worst possible score
// The instance should be feasible
// Get the cost (number of selected time-intervals) of a solution
// This method uses dynamic programming to solve the problem fast
int SolveFast(
    const int t,
    const int n,
    const int m,
    std::vector<Edge> &edges,
    bool (*could_be_specific_step_check)(int, std::vector<Edge>&, std::set<int>&, std::set<int>&, std::vector<int>&, std::vector<std::vector<Interval>::iterator>&),
    const bool compute_best) {
  // cur_intervals_its[i]: the iterator to the i-th edge's first time-interval that ends at or after time point cur_t.
  // If there is no such time-interval then cur_intervals_its[i] == edges[i].intervals.end()
  std::vector<std::vector<Interval>::iterator> cur_intervals_its;
  for (int edge_id = 0; edge_id < m; ++edge_id) {
    cur_intervals_its.push_back(edges[edge_id].intervals.begin());
  }
  std::vector<int> prev_available_edges_ids;  // The ids of the edges which have an interval that covers time point (cur_t - 1)
  // The min cost (number of selected edges) to cover all time points in [0, cur_t - 1]
  // where the selected edges of time point (cur_t - 1) are the combination specified by the index in the vector.
  // The i-th index is the i-th valid (covering) edge combination (binary counting)
  std::vector<int> prev_dp_min_costs;
  std::vector<std::vector<std::set<int>>> prev_covering_combinations;
  prev_dp_min_costs.push_back(0);
  prev_covering_combinations.push_back(std::vector<std::set<int>>(1));
  std::vector<unsigned int> argmin_prev_covering_combination_idxs;  // Vector of size t (first item will be 0). Can be used to reconstruct an optimal solution.
  for (int cur_t = 0; cur_t < t; ++cur_t) {
    std::vector<int> available_edges_ids;  // The ids of the edges which have an interval that covers time point cur_t
    for (int edge_id = 0; edge_id < m; ++edge_id) {
      if (cur_intervals_its[edge_id] != edges[edge_id].intervals.end() && cur_intervals_its[edge_id]->start <= cur_t) {
        available_edges_ids.push_back(edge_id);
      }
    }
    // The min cost (number of selected edges) to cover all time points in [0, cur_t]
    // where the selected edges of time point cur_t are the combination specified by the index in the vector.
    // The i-th index is the i-th valid (covering) edge combination (binary counting)
    std::vector<int> dp_min_costs;
    std::vector<std::set<int>> covering_combinations;
    SubsetIterator<int> selectedAvailableEdgesSubsetIterator(available_edges_ids);
    std::vector<std::set<int>> relevant_edges_ids_vec;
    for (std::set<int> prev_covering_combination : prev_covering_combinations.back()) {
      relevant_edges_ids_vec.push_back(GetRelevantEdgesIdsSubset(edges, prev_covering_combination, cur_intervals_its, cur_t));
    }
    while (!selectedAvailableEdgesSubsetIterator.HasReachedEndOfIterator()) {
      std::set<int> covering_combination = selectedAvailableEdgesSubsetIterator.GetSubset();
      ConnectedComponentsDFSController cc_controller(n);
      for (int edge_id : covering_combination) {
        cc_controller.Connect(edges[edge_id].u, edges[edge_id].v);
      }
      if (cc_controller.IsConnected()) {
        int dp_min_cost = -1;
        unsigned int argmin_prev_covering_combination_idx = 0;  // Initialize to suppress warning
        for (unsigned int prev_covering_combination_idx = 0; prev_covering_combination_idx < prev_dp_min_costs.size(); ++prev_covering_combination_idx) {
          // There is no need to check if the prev_covering_combination contains a time-interval that is not included in the covering_combination
          // because if this is the case then a prev_covering_combination with this edge will be at least as good
          std::set<int> relevant_edges_ids = relevant_edges_ids_vec[prev_covering_combination_idx];
          const bool could_be_general_step = CouldBeGeneralStepCheck(relevant_edges_ids, covering_combination);
          const bool could_be_specific_step = could_be_specific_step_check(n, edges, relevant_edges_ids, covering_combination, available_edges_ids, cur_intervals_its);
          if (could_be_general_step && could_be_specific_step) {
            const int num_new_edges = CountNewEdges(relevant_edges_ids, covering_combination);
            const int cost_using_prev_combination = prev_dp_min_costs[prev_covering_combination_idx] + num_new_edges;
            if (dp_min_cost == -1 || (compute_best && cost_using_prev_combination < dp_min_cost) || (!compute_best && cost_using_prev_combination > dp_min_cost)) {
              dp_min_cost = cost_using_prev_combination;
              argmin_prev_covering_combination_idx = prev_covering_combination_idx;
            }
          }
        }
        if (dp_min_cost != -1) {
          dp_min_costs.push_back(dp_min_cost);
          argmin_prev_covering_combination_idxs.push_back(argmin_prev_covering_combination_idx);
          covering_combinations.push_back(covering_combination);
        }
      }
    }
    prev_available_edges_ids = available_edges_ids;
    prev_dp_min_costs = dp_min_costs;
    prev_covering_combinations.push_back(covering_combinations);
    // Update cur_intervals_its
    for (int edge_id = 0; edge_id < m; ++edge_id) {
      if (cur_intervals_its[edge_id] != edges[edge_id].intervals.end() && cur_intervals_its[edge_id]->end < cur_t + 1) {
        ++cur_intervals_its[edge_id];
      }
    }
  }
  if (compute_best) {
    return *min_element(prev_dp_min_costs.begin(), prev_dp_min_costs.end());
  }
  return *max_element(prev_dp_min_costs.begin(), prev_dp_min_costs.end());
}
