add_executable(evolve main.cpp pauli.cpp hamiltonian.cpp calculator.cpp)

target_link_libraries(evolve PRIVATE ginac::ginac cxxopts)
target_compile_options(
  evolve PRIVATE $<$<CONFIG:Debug>: -ggdb3> $<$<CONFIG:Release>: -flto=auto
                 -march=native -mtune=native>)
