add_executable(foo main.cpp pauli.cpp hamiltonian.cpp calculator.cpp)

target_link_libraries(foo PRIVATE ginac::ginac)
