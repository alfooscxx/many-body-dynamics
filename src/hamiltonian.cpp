#include "hamiltonian.h"

#include <algorithm>
#include <numeric>

#include "pauli.h"

void hamiltonian::group::do_coloring() {
  pauli_string::qubit_mask_t mask{};
  for (const auto& [string, _] : base_strings_) {
    mask |= string.sites();
  }
  const auto sites = mask_to_vector(mask);
  starting_point = sites.front();
  if (sites.size() > 1) {
    block_size_ = sites[1] - starting_point;
    for (auto it = std::next(sites.begin()); it != sites.end(); ++it) {
      block_size_ = std::gcd(block_size_, *it - *std::prev(it));
      if (block_size_ == 1) {
        break;
      }
    }
  }
  period_length_ = sites.back() - starting_point;
  period_length_ /= block_size_;
  ++period_length_;
}

[[nodiscard]] pauli_string_combination hamiltonian::group::filter(
    std::size_t color, int site) const {  // NOLINT
  pauli_string_combination result;
  for (const auto& [string, coef] : base_strings_) {
    for (int string_site : mask_to_vector(string.sites())) {
      int shift = site - string_site;
      std::size_t shift_color = color_rule(shift);
      if (shift_color == color) {
        result.emplace(string.translate(shift), coef);
      }
    }
  }
  return result;
}

[[nodiscard]] std::size_t hamiltonian::group::color_rule(int shift) const {
  if (shift >= 0) {
    return (shift / block_size_) % period_length_;
  }
  shift = -shift - 1;
  return period_length_ - 1 - shift / block_size_;
}

void hamiltonian::group_by_commutativity(pauli_string_combination&& sum) {
  const size_t graph_size = sum.size();
  std::vector<std::vector<int>> adjacency_matrix(graph_size,
                                                 std::vector<int>(graph_size));
  {
    auto iit = sum.begin();
    for (size_t i = 0; i < graph_size; ++i, ++iit) {
      auto jit = sum.begin();
      for (size_t j = 0; j < graph_size; ++j, ++jit) {
        adjacency_matrix[i][j] =
            static_cast<int>(!iit->first.does_commute_with(jit->first));
      }
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
  std::vector<int> groups(graph_size, -1);
  auto can_merge_with = [&adjacency_matrix, &groups](size_t vertex,
                                                     size_t try_group) {
    for (size_t i = 0; i < adjacency_matrix.size(); ++i) {
      if (adjacency_matrix[vertex][i] != 0 && try_group == groups[i]) {
        return false;
      }
    }
    return true;
  };

  int chi = -1;
  for (size_t i = 0; i < graph_size; ++i) {
    size_t vertex = order[i];
    for (groups[vertex] = 0; !can_merge_with(vertex, groups[vertex]);
         ++groups[vertex]) {
      chi = std::max(chi, groups[vertex]);
    }
  }
  base_strings_groups_.resize(chi + 2);
  std::vector<pauli_string_combination> result(chi + 2);
  {
    auto iit = sum.begin();
    for (size_t i = 0; i < graph_size; ++i, ++iit) {
      base_strings_groups_[groups[order[i]]].emplace(iit->first, iit->second);
    }
  }
}
