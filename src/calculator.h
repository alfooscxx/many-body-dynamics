#pragma once

#include <ginac/symbol.h>

#include <ranges>

#include "hamiltonian.h"
#include "pauli.h"

class evolution_calculator {
 public:
  using state = std::vector<std::pair<pauli_string, GiNaC::ex>>;

  evolution_calculator(const scaled_pauli_string& observable,
                       hamiltonian&& hamiltonian);

  void advance(std::size_t count = 1);

  auto show() { return std::views::all(state_); }

  auto show_strings() { return std::views::keys(state_); }

  GiNaC::possymbol& get_tau() { return tau_; }

 private:
  GiNaC::possymbol tau_;
  hamiltonian hamiltonian_;
  state state_;
  state new_state_;
  size_t n_{};
};
