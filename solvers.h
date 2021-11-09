#ifndef MST_SOLVERS_H_
#define MST_SOLVERS_H_

#include <vector>
#include <set>

#include "general.h"

void InitRand();

// t: number of time points
// n: number of nodes
// m: number of edges
// edges: a vector with the graph's edges
// Check if an instance is feasible
bool IsFeasibleCheck(int t, int n, int m, std::vector<Edge> &edges);

// t: number of time points
// n: number of nodes
// m: number of edges
// edges: a vector with the graph's edges
// Get the cost (number of selected time-intervals) of the greedy algorithm's solution, or -1 if the instance is infeasible
int SolveGreedy(int t, int n, int m, std::vector<Edge> &edges);

// t: number of time points
// n: number of nodes
// m: number of edges
// edges: a vector with the graph's edges
// Get the cost (number of selected time-intervals) of an optimal solution, or -1 if the instance is infeasible
int SolveOptimal(int t, int n, int m, std::vector<Edge> &edges);

bool CouldBeOptimalStepCheck(
    int n,
    std::vector<Edge> &edges,
    std::set<int> &relevant_edges_ids,
    std::set<int> &covering_combination,
    std::vector<int> &available_edges_ids,
    std::vector<std::vector<Interval>::iterator> &cur_intervals_its);

bool CouldBeGreedyStepCheck(
    int n,
    std::vector<Edge> &edges,
    std::set<int> &relevant_edges_ids,
    std::set<int> &covering_combination,
    std::vector<int> &available_edges_ids,
    std::vector<std::vector<Interval>::iterator> &cur_intervals_its);

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
    int t,
    int n,
    int m,
    std::vector<Edge> &edges,
    bool (*could_be_specific_step_check)(int, std::vector<Edge>&, std::set<int>&, std::set<int>&, std::vector<int>&, std::vector<std::vector<Interval>::iterator>&),
    bool compute_best);

#endif  // MST_SOLVERS_H_
