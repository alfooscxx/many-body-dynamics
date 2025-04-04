#pragma once

#include <limits>
#include <span>

#include "pauli.h"

class hamiltonian {
 public:
  class group {
   public:
    [[nodiscard]] std::size_t color_number() const { return period_length_; }

    void do_coloring();

    void emplace(pauli_string string, const GiNaC::ex& coef) {
      base_strings_.emplace(string, coef);
    }

    [[nodiscard]] pauli_string_combination filter(std::size_t color,
                                                  int site) const;

   private:
    [[nodiscard]] std::size_t color_rule(int shift) const;

    pauli_string_combination base_strings_;
    // coloring data
    int starting_point{std::numeric_limits<int>::max()};
    std::size_t block_size_{1};
    std::size_t period_length_{1};
  };

  explicit hamiltonian(pauli_string_combination&& cycle_term) {
    group_by_commutativity(std::move(cycle_term));
    for (auto& group : base_strings_groups_) {
      group.do_coloring();
    }
  }

  [[nodiscard]] std::span<group, std::dynamic_extent> groups_view() {
    return {base_strings_groups_};
  }

 private:
  void group_by_commutativity(pauli_string_combination&&);

  std::vector<group> base_strings_groups_;
};
