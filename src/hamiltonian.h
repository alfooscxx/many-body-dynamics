#pragma once

#include "pauli.h"

class hamiltonian {
 public:
  explicit hamiltonian(pauli_string_combination&& sum) : sum(std::move(sum)) {}

  std::vector<pauli_string_combination> group_by_commutativity();

 private:
  pauli_string_combination sum;
};
