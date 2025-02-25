#include "calculator.h"

#include <string>

#include "inifcns.h"
#include "operators.h"
#include "pauli.h"
#include "relational.h"
#include "symbol.h"
#include "utils.h"

evolution_calculator::evolution_calculator(GiNaC::ex observable,
                                           const hamiltonian& hamiltonian)
    : state_(std::move(observable)), generic_evolution_operator_(GiNaC::_ex1) {
  const auto groups = hamiltonian.group_by_commutativity();
  for (const auto& H_i : groups) {
    for (const auto& term : H_i) {
      generic_evolution_operator_ *= time_evolution(term, tau);
    }
  }
}

void evolution_calculator::advance(size_t count) {
  for (size_t i = 0; i < count; ++i) {
    ++n_;
    GiNaC::possymbol tau_n("tau_" + std::to_string(n_));
    const GiNaC::ex evolution_operator =
        generic_evolution_operator_.subs(tau == tau_n);
    state_ = evolution_operator * state_ * evolution_operator.conjugate();
  }
}

GiNaC::ex evolution_calculator::time_evolution(scaled_pauli_string term,
                                               const GiNaC::ex& tau) {
  const auto [P, coef] = std::move(term);
  const auto phase_adjustment = P.phase_adjustment();
  const auto arg = phase_adjustment * coef * tau;
  return (-phase_adjustment) *
         (GiNaC::cos(arg) * pauli_string{} + GiNaC::I * GiNaC::sin(arg) * P);
}
