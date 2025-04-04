#pragma once

#include <initializer_list>

#include "basic.h"
#include "ex.h"
#include "flags.h"
#include "lst.h"
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
 public:
  using qubit_mask_t = unsigned long;
  enum class pauli_matrix : uint64_t { ONE = 0, X, Z, Y };
  friend scaled_pauli_string make_pauli_string(
      std::initializer_list<std::pair<int, pauli_matrix>>);
  friend bool operator<(const pauli_string& lhs, const pauli_string& rhs) {
    return std::tie(lhs.v_, lhs.w_) < std::tie(rhs.v_, rhs.w_);
  }

  // provides default ctor, duplicate(), accept(visitor&)
  // and compare_same_type(const basic&)
  GINAC_DECLARE_REGISTERED_CLASS(pauli_string, GiNaC::basic);  // NOLINT

  // for internal use only
  // creates pauli string with a single pauli matrix at the specified site
  // phase-wise
  explicit pauli_string(std::size_t site, pauli_matrix matrix);

 public:
  bool is_zero() const;

  // checks if two pauli strings commute with each other
  bool does_commute_with(const pauli_string& other) const;

  pauli_string translate(int shift) const;

  qubit_mask_t sites() const { return v_ | w_; };

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

  std::pair<qubit_mask_t, qubit_mask_t> representation() const {
    return {v_, w_};
  }

  qubit_mask_t v_{};
  qubit_mask_t w_{};
};

std::vector<int> mask_to_vector(pauli_string::qubit_mask_t);

struct scaled_pauli_string {
  pauli_string P;
  GiNaC::ex coef;
};

using pauli_string_combination = std::map<pauli_string, GiNaC::ex>;

class find_all_pauli_strings : public GiNaC::visitor,
                               public pauli_string::visitor {
  GiNaC::lst list_;
  void visit(const pauli_string& string) override { list_.append(string); }

 public:
  const GiNaC::lst& get_result() {
    list_.sort();
    list_.unique();
    return list_;
  }
};
