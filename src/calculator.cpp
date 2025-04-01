#include "calculator.h"

#include "hamiltonian.h"
#include "pauli.h"

using namespace GiNaC;

evolution_calculator::evolution_calculator(
    const scaled_pauli_string& observable, hamiltonian&& hamiltonian)
    : hamiltonian_(std::move(hamiltonian)) {
  state_.emplace_back(observable.P, observable.coef);
}

void evolution_calculator::advance(std::size_t count) {
  static const auto arg_coef = 2 * tau_;
  for (std::size_t iter = 0; iter < count; ++iter) {
    ++n_;
    for (const auto& group : hamiltonian_.groups_view()) {
      for (std::size_t color = 0; color < group.color_number(); ++color) {
        pauli_string::qubit_mask_t mask{};
        for (const auto& string : std::views::keys(state_)) {
          mask |= string.sites();
        }
        const auto sites = mask_to_vector(mask);
        pauli_string_combination conflicts;
        for (int site : sites) {
          for (const auto& [string, coef] : group.filter(color, site)) {
            conflicts.try_emplace(string, coef);
          }
        }
        for (const auto& [P, P_coef] : conflicts) {
          for (const auto& [A, A_coef] : state_) {
            if (P.does_commute_with(A)) {
              new_state_.emplace_back(A, A_coef);
              continue;
            }
            // cos(cP) = cos(c)
            // sin(cP) = sin(c) * P
            const auto P_phase_adjustment = P.phase_adjustment();
            const auto arg = arg_coef * P_phase_adjustment * P_coef;
            const auto [PA, PA_phase_adjustment] = P * A;
            new_state_.emplace_back(A, GiNaC::cos(arg) * A_coef);
            new_state_.emplace_back(
                PA, PA_phase_adjustment * P_phase_adjustment.conjugate() *
                        GiNaC::I * GiNaC::sin(arg) * A_coef);
          }
          std::ranges::sort(new_state_, {},
                            &decltype(new_state_)::value_type::first);
          for (std::size_t i = new_state_.size() - 1; i > 0; --i) {
            if (new_state_[i].first == new_state_[i - 1].first) {
              new_state_[i - 1].second += new_state_[i].second;
            }
          }
          const auto junk = std::ranges::unique(
              new_state_, {}, &decltype(new_state_)::value_type::first);
          new_state_.erase(junk.begin(), junk.end());
          const auto removed = std::ranges::remove_if(
              new_state_, [](const auto& coef) { return coef.is_zero(); },
              &decltype(new_state_)::value_type::second);
          new_state_.erase(removed.begin(), removed.end());
          state_.swap(new_state_);
          new_state_.clear();
        }
      }
    }
  }
}
