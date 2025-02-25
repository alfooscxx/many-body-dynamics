#pragma once

#include "basic.h"
#include "ex.h"
#include "flags.h"
#include "print.h"
#include "registrar.h"

struct scaled_pauli_string;

/**
 * @brief
 * Pauli string is a multi-qubit operator
 * represented by a tensor product of several pauli matrices acting on different
 * sites
 */
class pauli_string : public GiNaC::basic {
  using inner_bitstring = unsigned long;

 public:
  enum class pauli_matrix : uint64_t { ONE = 0, X, Z, Y };
  friend scaled_pauli_string make_pauli_string(std::size_t site,
                                               pauli_matrix matrix);

  // provides default ctor, duplicate(), accept(visitor&)
  // and compare_same_type(const basic&)
  GINAC_DECLARE_REGISTERED_CLASS(pauli_string, GiNaC::basic);  // NOLINT

  // for internal use only
  // creates pauli string with a single pauli matrix at the specified site
  // phase-wise
  explicit pauli_string(std::size_t site, pauli_matrix matrix);

 public:
  // checks if two pauli strings commute with each other
  bool does_commute_with(const pauli_string& other) const;

  GiNaC::ex phase_adjustment() const;

 protected:
  bool is_equal_same_type(const GiNaC::basic& other) const override;
  GiNaC::ex eval_ncmul(const GiNaC::exvector& mul) const override;
  unsigned return_type() const override {
    return GiNaC::return_types::noncommutative;
  }
  GiNaC::return_type_t return_type_tinfo() const override;
  unsigned int calchash() const override;
  void do_print(const GiNaC::print_context& context, unsigned int level) const;

 private:
  static unsigned int get_hash_seed();

  std::pair<inner_bitstring, inner_bitstring> representation() const {
    return {v_, w_};
  }

  inner_bitstring v_{};
  inner_bitstring w_{};
};

struct scaled_pauli_string {
  pauli_string P;
  GiNaC::ex coef;
};

using pauli_string_combination = std::vector<scaled_pauli_string>;
