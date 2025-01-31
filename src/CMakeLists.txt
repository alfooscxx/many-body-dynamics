add_executable(foo main.cpp pauli.cpp)

target_link_libraries(foo PRIVATE ginac::ginac)
