#pragma once

#include <ginac/ginac.h>

struct scaled_pauli_string;

/**
 * @brief
 * Pauli string is a multi-qubit operator
 * represented by a tensor product of several pauli matrices acting on different
 * sites
 */
class pauli_string {
 public:
  using qubit_mask_t = unsigned long;
  enum class pauli_matrix : uint8_t { ONE = 0, X, Z, Y };

  friend scaled_pauli_string make_pauli_string(
      const std::vector<std::pair<int, pauli_matrix>>&);
  friend auto operator<=>(const pauli_string& lhs,
                          const pauli_string& rhs) = default;
  friend scaled_pauli_string operator*(const pauli_string& lhs,
                                       const pauli_string& rhs);
  friend std::ostream& operator<<(std::ostream& stream,
                                  const pauli_string& string);

  pauli_string() = default;

  // for internal use only
  // creates pauli string with a single pauli matrix at the specified site
  // phase-wise
  explicit pauli_string(std::size_t site, pauli_matrix matrix);

  // checks if two pauli strings commute with each other
  [[nodiscard]] bool does_commute_with(const pauli_string& other) const;

  [[nodiscard]] pauli_string translate(int shift) const;

  [[nodiscard]] qubit_mask_t sites() const { return v_ | w_; };

  [[nodiscard]] GiNaC::numeric phase_adjustment() const;

  [[nodiscard]] GiNaC::numeric polarize(GiNaC::numeric p_x, GiNaC::numeric p_y,
                                        GiNaC::numeric p_z) const;

 private:
  [[nodiscard]] std::pair<qubit_mask_t, qubit_mask_t> representation() const {
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
