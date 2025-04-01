#include "pauli.h"

#include <bit>
#include <climits>
#include <utility>

scaled_pauli_string make_pauli_string(
    const std::vector<std::pair<int, pauli_string::pauli_matrix>>& init) {
  pauli_string string;
  for (auto [site, matrix] : init) {
    pauli_string temp{static_cast<size_t>(site), matrix};
    string.v_ ^= temp.v_;
    string.w_ ^= temp.w_;
  }
  return {.P = string, .coef = string.phase_adjustment().conjugate()};
}

pauli_string::pauli_string(std::size_t site, pauli_string::pauli_matrix matrix)
    : pauli_string() {
  const auto matrix_byte = std::to_underlying(matrix);
  // clang-format off
  v_ = ((matrix_byte & 2ULL) >> 1ULL) << site;
  w_ = ( matrix_byte & 1ULL)          << site;
  // clang-format on
}

namespace {
template <typename T>
unsigned int popcount(T value) {
  if constexpr (std::is_integral_v<T>) {
    return std::popcount(value);
  } else {
    return value.count();
  }
}
}  // namespace

scaled_pauli_string operator*(const pauli_string& lhs,
                              const pauli_string& rhs) {
  pauli_string result;
  result.v_ = lhs.v_ ^ rhs.v_;
  result.w_ = lhs.w_ ^ rhs.w_;
  return {.P = result,
          .coef = (((popcount(lhs.w_ & rhs.v_) % 2) == 0U) ? 1 : -1)};
}

bool pauli_string::does_commute_with(const pauli_string& other) const {
  return (popcount(v_ & other.w_) ^ popcount(w_ & other.v_)) % 2 == 0;
}

pauli_string pauli_string::translate(int shift) const {
  pauli_string result;
  if (shift >= 0) {
    result.v_ = v_ << static_cast<unsigned>(shift);
    result.w_ = w_ << static_cast<unsigned>(shift);
  } else {
    result.v_ = v_ >> static_cast<unsigned>(-shift);
    result.w_ = w_ >> static_cast<unsigned>(-shift);
  }
  return result;
}

std::vector<int> mask_to_vector(pauli_string::qubit_mask_t sites) {
  std::vector<int> result;
  while (sites != 0) {
    result.push_back(std::countr_zero(sites));
    sites &= (sites - 1);
  }
  return result;
}

GiNaC::numeric pauli_string::phase_adjustment() const {
  static const std::array<GiNaC::numeric, 4> phase = {1, GiNaC::I, -1,
                                                      -GiNaC::I};
  return phase[popcount(v_ & w_) % 4];
}

[[nodiscard]] GiNaC::numeric pauli_string::polarize(GiNaC::numeric p_x,
                                                    GiNaC::numeric p_y,
                                                    GiNaC::numeric p_z) const {
  GiNaC::numeric result = phase_adjustment();
  std::array<GiNaC::numeric, 4> substitution = {1, std::move(p_x),
                                                std::move(p_z), std::move(p_y)};
  const std::size_t bits = sizeof(qubit_mask_t) * CHAR_BIT;
  for (std::size_t i = 0; i < bits; ++i) {
    const auto matrix_byte =
        (((v_ & (1ULL << i)) >> i) << 1ULL) | ((w_ & (1ULL << i)) >> i);
    result *= substitution.at(matrix_byte);
  }
  return result;
}

std::ostream& operator<<(std::ostream& stream, const pauli_string& string) {
  auto [v_, w_] = string.representation();
  static constexpr std::array print_matrix = {"", "X", "Z", "Y"};
  static constexpr std::array print_phase = {"", "I", "-", "-I"};
  if ((v_ == 0U) && (w_ == 0U)) {
    stream << "ONE";
    return stream;
  }
  const std::size_t bits = sizeof(pauli_string::qubit_mask_t) * CHAR_BIT;
  const auto phase_degree = popcount(v_ & w_) % 4;
  stream << print_phase[phase_degree];  // NOLINT
  for (std::size_t i = 0; i < bits; ++i) {
    const auto matrix_byte =
        (((v_ & (1ULL << i)) >> i) << 1ULL) | ((w_ & (1ULL << i)) >> i);
    if (matrix_byte != 0) {
      stream << "[" << print_matrix.at(matrix_byte) << "_" << i << "]";
    }
  }
  return stream;
}
