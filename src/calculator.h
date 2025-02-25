#pragma once

#include "hamiltonian.h"
#include "pauli.h"
#include "symbol.h"

class evolution_calculator {
 public:
  evolution_calculator(GiNaC::ex observable, const hamiltonian& hamiltonian);

  void advance(size_t count = 1);

  GiNaC::ex show() { return state_; }

 private:
  static GiNaC::ex time_evolution(scaled_pauli_string term,
                                  const GiNaC::ex& tau);

  GiNaC::ex state_;
  GiNaC::ex generic_evolution_operator_;
  GiNaC::possymbol tau;
  size_t n_{};
};
