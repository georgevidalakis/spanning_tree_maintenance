#include <iostream>
#include <vector>

#include "general.h"
#include "graph_sampler.h"
#include "solvers.h"

int main() {
  InitRand();
  int num_trials;
  std::cin >> num_trials;
  int t, d, n, m;
  std::cin >> t >> d >> n >> m;
  GraphSampler graph_sampler(t, d, n, m);
  for (int edge_id = 0; edge_id < m; ++edge_id) {
    int u, v;
    std::cin >> u >> v;
    graph_sampler.AddEdge(u, v);
  }
  int num_infeasible = 0, num_greedy_is_optimal = 0, num_exists_optimal_greedy = 0, num_greedy_is_suboptimal = 0;
  int max_ratio_nom = 0, max_ratio_denom = 1;
  for (int trial_id = 0; trial_id < num_trials; ++trial_id) {
    std::vector<Edge> edges = graph_sampler.GetEdgesWithRefreshedIntervals();
    bool is_feasible = IsFeasibleCheck(t, n, m, edges);
    if (!is_feasible) {
      ++num_infeasible;
      continue;
    }
    const int fast_optimal_cost = SolveFast(t, n, m, edges, CouldBeOptimalStepCheck, true);
    const int fast_best_greedy_cost = SolveFast(t, n, m, edges, CouldBeGreedyStepCheck, true);
    const int fast_worst_greedy_cost = SolveFast(t, n, m, edges, CouldBeGreedyStepCheck, false);
    if (fast_optimal_cost > fast_best_greedy_cost || fast_best_greedy_cost > fast_worst_greedy_cost) {
      std::cerr << "[BUG]:\n";
      std::cerr << "fast_optimal_cost: " << fast_optimal_cost << '\n';
      std::cerr << "fast_best_greedy_cost: " << fast_best_greedy_cost << '\n';
      std::cerr << "fast_worst_greedy_cost: " << fast_worst_greedy_cost << '\n';
      exit(1);
    }
    /*// Sanity check (slow)
    const int optimal_cost = SolveOptimal(t, n, m, edges);
    const int greedy_cost = SolveGreedy(t, n, m, edges);
    if (optimal_cost != fast_optimal_cost || fast_best_greedy_cost > greedy_cost || greedy_cost > fast_worst_greedy_cost) {
      std::cerr << "[BUG]:\n";
      std::cerr << "greedy_cost: " << greedy_cost << '\n';
      std::cerr << "optimal_cost: " << optimal_cost << '\n';
      std::cerr << "fast_optimal_cost: " << fast_optimal_cost << '\n';
      std::cerr << "fast_best_greedy_cost: " << fast_best_greedy_cost << '\n';
      std::cerr << "fast_worst_greedy_cost: " << fast_worst_greedy_cost << '\n';
      exit(1);
    }*/
    if (fast_worst_greedy_cost == fast_optimal_cost) {
      ++num_greedy_is_optimal;
    }
    else if (fast_best_greedy_cost == fast_optimal_cost) {
      ++num_exists_optimal_greedy;
    }
    else {
      ++num_greedy_is_suboptimal;
      /*std::cout << "SUBOPT:\n";
      std::cout << "On trial " << trial_id << ":\n";
      for (int edge_id = 0; edge_id < m; ++edge_id) {
        std::cout << "Edge " << edge_id << ":\n";
        for (Interval interval : edges[edge_id].intervals) {
          std::cout << '[' << interval.start << ", " << interval.end << "]\n";
        }
        std::cout << '\n';
      }
      std::cout << '\n';*/
    }
    if (max_ratio_nom * fast_optimal_cost < max_ratio_denom * fast_worst_greedy_cost) {
      max_ratio_nom = fast_worst_greedy_cost;
      max_ratio_denom = fast_optimal_cost;
    }
  }
  // Output experiments' statistics
  std::cout << "Num infeasible: " << num_infeasible << '\n';
  std::cout << "Num greedy is OPT: " << num_greedy_is_optimal << '\n';
  std::cout << "Num exists OPT greedy: " << num_exists_optimal_greedy << '\n';
  std::cout << "Num greedy is SUBOPT: " << num_greedy_is_suboptimal << '\n';
  std::cout << '\n';
  // Output max ratio
  std::cout << "Max ratio: " << max_ratio_nom << '/' << max_ratio_denom << '\n';
  return 0;
}
