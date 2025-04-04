#pragma once

#include "hamiltonian.h"
#include "pauli.h"
#include "symbol.h"

class evolution_calculator {
 public:
  evolution_calculator(GiNaC::ex observable, hamiltonian&& hamiltonian);

  void advance(size_t count = 1);

  GiNaC::ex show() { return state_; }

 private:
  static GiNaC::ex time_evolution(scaled_pauli_string term,
                                  const GiNaC::ex& tau);

  static const GiNaC::exmap simplify_exponents;

  hamiltonian hamiltonian_;
  GiNaC::ex state_;
  GiNaC::possymbol tau_{"tau"};
  size_t n_{};
};
