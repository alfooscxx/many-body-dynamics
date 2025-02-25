#include "hamiltonian.h"

#include <algorithm>
#include <numeric>

#include "pauli.h"

std::vector<pauli_string_combination> hamiltonian::group_by_commutativity()
    const {
  const size_t graph_size = sum.size();
  std::vector<std::vector<int>> adjacency_matrix(graph_size,
                                                 std::vector<int>(graph_size));
  for (size_t i = 0; i < graph_size; ++i) {
    for (size_t j = 0; j < graph_size; ++j) {
      adjacency_matrix[i][j] =
          static_cast<int>(!sum[i].P.does_commute_with(sum[j].P));
    }
  }
  std::vector<unsigned long long> order(graph_size);
  std::ranges::iota(order, 0);
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
    for (size_t i = 0; i < adjacency_matrix.size(); ++i) {
      if (adjacency_matrix[vertex][i] != 0 && try_color == color[i]) {
        return false;
      }
    }
    return true;
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
