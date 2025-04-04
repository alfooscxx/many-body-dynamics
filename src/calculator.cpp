#include "calculator.h"

#include "ex.h"
#include "ginac.h"
#include "hamiltonian.h"
#include "inifcns.h"
#include "operators.h"
#include "pauli.h"
#include "relational.h"
#include "symbol.h"
#include "utils.h"
#include "wildcard.h"

using namespace GiNaC;

ex simplify(const ex& expr, const exmap& rules) {
  if (is_a<mul>(expr)) {
    ex next = expr;
    ex prev = next;
    next = next.subs(rules, subs_options::subs_algebraic).expand();
    while (next != prev) {
      prev = next;
      next = next.subs(rules, subs_options::subs_algebraic).expand();
    }
    return next;
  }
  if (is_a<add>(expr)) {
    exvector collected;
    for (size_t i = 0; i < expr.nops(); ++i) {
      collected.push_back(simplify(expr.op(i), rules));
    }
    return add(collected);
  }
  return expr;
}

// NOLINTBEGIN
const GiNaC::exmap evolution_calculator::simplify_exponents{
    {exp(wild(0)) * exp(wild(1)), exp(wild(0) + wild(1))}};
// NOLINTEND

evolution_calculator::evolution_calculator(GiNaC::ex observable,
                                           hamiltonian&& hamiltonian)
    : state_(std::move(observable)), hamiltonian_(std::move(hamiltonian)) {}

void evolution_calculator::advance(size_t count) {
  for (size_t iter = 0; iter < count; ++iter) {
    ++n_;
    // GiNaC::possymbol tau_n("tau_" + std::to_string(n_));
    for (const auto& group : hamiltonian_.groups_view()) {
      for (std::size_t color = 0; color < group.color_number(); ++color) {
        find_all_pauli_strings visitor;
        state_.traverse_postorder(visitor);
        const auto strings_in_A = visitor.get_result();
        pauli_string::qubit_mask_t mask{};
        for (const auto& string : strings_in_A) {
          mask |= GiNaC::ex_to<pauli_string>(string).sites();
        }
        const auto sites = mask_to_vector(mask);
        pauli_string_combination conflicts;
        for (int site : sites) {
          for (const auto& [string, coef] : group.filter(color, site)) {
            bool non_commute = false;
            for (const auto& omega_string : strings_in_A) {
              if (!string.does_commute_with(
                      GiNaC::ex_to<pauli_string>(omega_string))) {
                non_commute = true;
                break;
              }
            }
            if (non_commute) {
              conflicts.emplace(string, coef);
            }
          }
        }
        for (const auto& [string, coef] : conflicts) {
          auto step = time_evolution({.P = string, .coef = coef}, tau_);
          state_ = step * state_ * step.conjugate();
          state_ = simplify(state_.expand(), simplify_exponents);
        }
      }
    }
  }
}

GiNaC::ex evolution_calculator::time_evolution(scaled_pauli_string term,
                                               const GiNaC::ex& tau) {
  const auto [P, coef] = std::move(term);
  const auto phase_adjustment = P.phase_adjustment();
  const auto arg = GiNaC::I * phase_adjustment * coef * tau;
  return (-phase_adjustment) *
         ((GiNaC::exp(arg) + GiNaC::exp(-arg)) / numeric(2) * pauli_string{} +
          (GiNaC::exp(arg) - GiNaC::exp(-arg)) / numeric(2) * P);
}
