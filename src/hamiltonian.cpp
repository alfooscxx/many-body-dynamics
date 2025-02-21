#include "hamiltonian.h"

#include <algorithm>
#include <numeric>
#include <ranges>

#include "pauli.h"

std::vector<pauli_string_combination> hamiltonian::group_by_commutativity() {
  const size_t graph_size = sum.size();
  std::vector<std::vector<int>> adjacency_matrix(graph_size,
                                                 std::vector<int>(graph_size));
  for (size_t i = 0; i < graph_size; ++i) {
    for (size_t j = 0; j < graph_size; ++j) {
      adjacency_matrix[i][j] =
          static_cast<int>(!sum[i].P.does_commute_with(sum[j].P));
    }
  }
  std::vector<unsigned long long> order =
      std ::views::iota(0ULL, graph_size) | std::ranges::to<std::vector>();
  std::vector<size_t> degrees(graph_size);
  std::ranges::transform(order, degrees.begin(), [&adjacency_matrix](int idx) {
    return std::accumulate(adjacency_matrix[idx].begin(),
                           adjacency_matrix[idx].end(), 0ULL);
  });
  std::ranges::sort(order, [&degrees](int first, int second) {
    return degrees[first] > degrees[second];
  });
  std::vector<int> color(graph_size, -1);
  auto can_color = [&adjacency_matrix, &color](size_t vertex,
                                               size_t try_color) {
    return std::ranges::all_of(adjacency_matrix[vertex],
                               [try_color](size_t neighbour_color) {
                                 return neighbour_color != try_color;
                               });
  };
  int chi = 0;
  for (size_t i = 0; i < graph_size; ++i) {
    size_t vertex = order[i];
    for (color[vertex] = 0; !can_color(vertex, color[vertex]);
         ++color[vertex]) {
      chi = std::max(chi, color[vertex]);
    }
  }
  std::vector<pauli_string_combination> result(chi + 1);
  for (size_t i = 0; i < graph_size; ++i) {
    result[color[order[i]]].push_back(sum[order[i]]);
  }
  return result;
}
