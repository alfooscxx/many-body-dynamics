add_executable(evolve main.cpp pauli.cpp hamiltonian.cpp calculator.cpp)

target_include_directories(evolve PRIVATE ${SYMENGINE_INCLUDE_DIRS})
target_link_libraries(evolve PRIVATE OpenMP::OpenMP_CXX ${SYMENGINE_LIBRARIES}
                                     cxxopts)
target_compile_options(
  evolve PRIVATE $<$<CONFIG:Debug>: -ggdb3> $<$<CONFIG:Release>: -flto=auto
                 -march=native -mtune=native>)
