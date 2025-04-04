#include "pauli.h"

#include <bit>
#include <climits>
#include <utility>

#include "hash_seed.h"
#include "ncmul.h"
#include "numeric.h"
#include "operators.h"
#include "utils.h"

scaled_pauli_string make_pauli_string(
    std::initializer_list<std::pair<int, pauli_string::pauli_matrix>> init) {
  pauli_string string;
  GiNaC::ex coef{1};
  for (auto [site, matrix] : init) {
    if (matrix == pauli_string::pauli_matrix::Y) {
      coef *= -GiNaC::I;
    }
    pauli_string temp{static_cast<size_t>(site), matrix};
    string.v_ ^= temp.v_;
    string.w_ ^= temp.w_;
  }
  return {.P = string, .coef = coef};
}

// NOLINTBEGIN
GINAC_IMPLEMENT_REGISTERED_CLASS_OPT(
    pauli_string, GiNaC::basic,
    print_func<GiNaC::print_context>(&pauli_string::do_print));
// NOLINTEND

pauli_string::pauli_string() {
  setflag(GiNaC::status_flags::evaluated | GiNaC::status_flags::expanded);
}

pauli_string::pauli_string(std::size_t site, pauli_string::pauli_matrix matrix)
    : pauli_string() {
  const auto matrix_byte = std::to_underlying(matrix);
  // clang-format off
  v_ = ((matrix_byte & 2ULL) >> 1ULL) << site;
  w_ = ( matrix_byte & 1ULL)          << site;
  // clang-format on
}

int pauli_string::compare_same_type(const GiNaC::basic& other) const {
  GINAC_ASSERT(is_a<pauli_string>(other));
  const auto& second = dynamic_cast<const pauli_string&>(other);
  if (v_ != second.v_) {
    return v_ < second.v_ ? -1 : 1;
  }
  if (w_ != second.w_) {
    return w_ < second.w_ ? -1 : 1;
  }
  return 0;
}

bool pauli_string::is_equal_same_type(const GiNaC::basic& other) const {
  GINAC_ASSERT(is_a<pauli_string>(other));
  const auto& second = dynamic_cast<const pauli_string&>(other);
  return std::tie(v_, w_) == std::tie(second.v_, second.w_);
}

template <typename T>
unsigned int popcount(T value) {
  if constexpr (std::is_integral_v<T>) {
    return std::popcount(value);
  } else {
    return value.count();
  }
}

bool pauli_string::is_zero() const { return (v_ | w_) == 0; }

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
  const std::size_t bits = sizeof(pauli_string::qubit_mask_t) * CHAR_BIT;
  std::vector<int> result;
  while (sites != 0) {
    result.push_back(std::countr_zero(sites));
    sites &= (sites - 1);
  }
  return result;
}

GiNaC::ex pauli_string::phase_adjustment() const {
  static const std::array<GiNaC::ex, 4> phase = {GiNaC::_ex1, GiNaC::I,
                                                 -GiNaC::_ex1, -GiNaC::I};
  return phase[popcount(v_ & w_) % 4];  // NOLINT
}

GiNaC::ex pauli_string::eval_ncmul(const GiNaC::exvector& mul) const {
  GiNaC::exvector result;
  int phase = 1;
  pauli_string current{};
  for (const auto& expr : mul) {
    if (!is_a<pauli_string>(expr)) {
      if ((current.v_ != 0U) || (current.w_ != 0U)) {
        result.emplace_back(current);
        current = pauli_string{};
      }
      result.push_back(expr);
      continue;
    }
    const auto [update_v, update_w] =
        GiNaC::ex_to<pauli_string>(expr).representation();
    if (popcount(current.v_ & update_w) % 2 == 1) {
      phase = -phase;
    }
    current.v_ ^= update_v;
    current.w_ ^= update_w;
  }
  if ((current.v_ != 0U) || (current.w_ != 0U)) {
    result.emplace_back(current);
  }
  if (result.empty()) {
    return phase * pauli_string{};
  }
  if (result.size() < mul.size()) {
    return phase * GiNaC::reeval_ncmul(result);
  }
  return phase * GiNaC::hold_ncmul(result);
}

GiNaC::return_type_t pauli_string::return_type_tinfo() const {
  return GiNaC::make_return_type_t<pauli_string>();
}

unsigned int pauli_string::get_hash_seed() {
  static const unsigned int hash_seed =
      GiNaC::make_hash_seed(typeid(pauli_string));
  return hash_seed;
}

unsigned int pauli_string::calchash() const {
  setflag(GiNaC::status_flags::hash_calculated);
  hashvalue = get_hash_seed();
  const auto V = static_cast<unsigned>(v_);  // NOLINT
  const auto W = static_cast<unsigned>(w_);  // NOLINT
  // an "elegant" integer pair hash function
  // http://szudzik.com/ElegantPairing.pdf
  if (V < W) {
    hashvalue ^= W * W + V;
    hashvalue = GiNaC::golden_ratio_hash(hashvalue);
    return hashvalue;
  }
  hashvalue ^= V * V + V + W;
  hashvalue = GiNaC::golden_ratio_hash(hashvalue);
  return hashvalue;
}

void pauli_string::do_print(const GiNaC::print_context& context,
                            unsigned int level) const {
  (void)level;
  static constexpr std::array print_matrix = {"", "X", "Z", "Y"};
  static constexpr std::array print_phase = {"", "I", "-", "-I"};
  if ((v_ == 0U) && (w_ == 0U)) {
    context.s << "ONE";
    return;
  }
  const std::size_t bits = sizeof(qubit_mask_t) * CHAR_BIT;
  const auto phase_degree = popcount(v_ & w_) % 4;
  context.s << print_phase[phase_degree];  // NOLINT
  for (std::size_t i = 0; i < bits; ++i) {
    const auto matrix_byte =
        (((v_ & (1ULL << i)) >> i) << 1ULL) | ((w_ & (1ULL << i)) >> i);
    if (matrix_byte != 0) {
      // NOLINTBEGIN
      context.s << "[" << print_matrix[matrix_byte] << "_" << i << "]";
      // NOLINTEND
    }
  }
}
