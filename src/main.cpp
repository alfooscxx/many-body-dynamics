// main.cpp — CLI tool for Suzuki–Trotter evolution of a quantum observable
#include <algorithm>
#include <cxxopts.hpp>
#include <iostream>
#include <numeric>
#include <ranges>
#include <sstream>
#include <string>
#include <vector>

#include "calculator.h"
#include "hamiltonian.h"
#include "pauli.h"

namespace rv = std::views;
namespace rg = std::ranges;

//--------------------------------------------------------------------
//   Command‑line options
//--------------------------------------------------------------------
struct Cli {
  int trotter_steps = 1;
  double dt = 0.1;     // time grid step
  double t_max = 1.0;  // final time
  std::array<double, 3> pol = {1.0, 0.0, 0.0};
  std::string hamiltonian = "XX+Z";
  std::string observable = "Z";
};

namespace {
Cli parse_cli(int argc, char** argv) {
  cxxopts::Options cli("trotter", "Suzuki–Trotter evolution simulator");
  cli.add_options()("steps", "Trotter steps",
                    cxxopts::value<int>()->default_value("1"))(
      "density", "Time grid step",
      cxxopts::value<double>()->default_value("0.1"))(
      "interval", "End of time interval",
      cxxopts::value<double>()->default_value("1.0"))(
      "substitution", "Polarization x,y,z",
      cxxopts::value<std::string>()->default_value("1,0,0"))(
      "hamiltonian", "Hamiltonian like XX+Z or XX+Z+X",
      cxxopts::value<std::string>()->default_value("XX+Z"))(
      "observable", "Observable like Z or XY",
      cxxopts::value<std::string>()->default_value("Z"))("help", "Show help");

  auto res = cli.parse(argc, argv);
  if (res.count("help")) {
    std::cout << cli.help() << '\n';
    std::exit(EXIT_SUCCESS);
  }

  Cli opt;
  opt.trotter_steps = res["steps"].as<int>();
  opt.dt = res["density"].as<double>();
  opt.t_max = res["interval"].as<double>();
  opt.hamiltonian = res["hamiltonian"].as<std::string>();
  opt.observable = res["observable"].as<std::string>();

  if (opt.trotter_steps < 1) throw std::invalid_argument("steps > 0 required");
  if (opt.dt <= 0) throw std::invalid_argument("density > 0 required");
  if (opt.t_max <= 0) throw std::invalid_argument("interval > 0 required");

  // parse polarization vector
  const auto csv = res["substitution"].as<std::string>();
  std::vector<double> vec;
  std::stringstream ss(csv);
  std::string token;
  while (std::getline(ss, token, ','))
    if (!token.empty()) vec.push_back(std::stod(token));
  if (vec.size() != 3)
    throw std::invalid_argument("substitution expects x,y,z");
  const double norm2 = vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2];
  if (std::abs(norm2 - 1.0) > 1e-6)
    throw std::invalid_argument("x^2+y^2+z^2 must equal 1");
  rg::copy(vec, opt.pol.begin());

  return opt;
}

scaled_pauli_string pauli_literal(std::string_view lit) {
  using M = pauli_string::pauli_matrix;
  std::vector<std::pair<int, M>> tmp;
  tmp.reserve(lit.size());
  for (auto [i, c] : rv::enumerate(lit)) {
    M m;
    switch (c) {
      case 'X':
        m = M::X;
        break;
      case 'Y':
        m = M::Y;
        break;
      case 'Z':
        m = M::Z;
        break;
      default:
        throw std::invalid_argument("invalid Pauli char");
    }
    tmp.emplace_back(static_cast<int>(i), m);
  }
  return make_pauli_string(tmp);
}
}  // namespace

//--------------------------------------------------------------------
//   Build Pauli string from literal like "XYZ"
//--------------------------------------------------------------------
int main(int argc, char** argv) {
  try {
    const Cli cfg = parse_cli(argc, argv);

    // Hamiltonian construction with ranges
    pauli_string_combination H_terms;
    auto add_term = [&](std::string_view part) {
      auto ps = pauli_literal(part);
      H_terms[ps.P] += ps.coef;  // default coefficient 1
    };
    rg::for_each(cfg.hamiltonian | rv::split('+'), [&](auto sub) {
      add_term({&*sub.begin(), static_cast<size_t>(rg::distance(sub))});
    });
    hamiltonian H(H_terms);

    // Observable
    auto observable = pauli_literal(cfg.observable);
    observable.P = observable.P.translate(32);

    // Evolution
    evolution_calculator calc(observable, std::move(H));
    calc.advance(cfg.trotter_steps);

    const auto& state = calc.show();

    for (double t = 0.0; t <= cfg.t_max + 1e-12; t += cfg.dt) {
      const SymEngine::Expression tau(t / cfg.trotter_steps);
      const SymEngine::Expression total = std::accumulate(
          state.begin(), state.end(), SymEngine::Expression{0},
          [&](const SymEngine::Expression& acc, const auto& term) {
            const auto& [P, coef] = term;
            const auto pol = P.polarize(cfg.pol[0], cfg.pol[1], cfg.pol[2]);
            return pol == 0 ? acc
                            : acc + coef.subs({{calc.get_tau(), tau}}) * pol;
          });
      std::cout << t << ' ' << static_cast<std::complex<double>>(total).real()
                << '\n';
    }
  } catch (const std::exception& e) {
    std::cerr << "Error: " << e.what() << '\n';
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
