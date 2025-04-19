find_package(OpenMP REQUIRED)

include(FetchContent)

FetchContent_Declare(
  symengine_src
  GIT_REPOSITORY https://github.com/symengine/symengine.git
  GIT_TAG v0.14.0)
FetchContent_Populate(symengine_src)

set(SYMENGINE_STAGE ${CMAKE_BINARY_DIR}/_symengine) # <build>/_symengine
set(SYMENGINE_BUILD ${CMAKE_BINARY_DIR}/_symengine-build) # separate build dir

if(NOT EXISTS ${SYMENGINE_BUILD}/build.ninja AND NOT EXISTS
                                                 ${SYMENGINE_BUILD}/Makefile)
  file(MAKE_DIRECTORY ${SYMENGINE_BUILD})
  execute_process(
    COMMAND
      ${CMAKE_COMMAND} -S ${symengine_src_SOURCE_DIR} -B ${SYMENGINE_BUILD}
      -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=${SYMENGINE_STAGE}
      -DBUILD_SHARED_LIBS=ON -DBUILD_TESTS=OFF -DBUILD_BENCHMARKS=OFF
      -DWITH_TCMALLOC=ON -DWITH_MPFR=ON -DWITH_LLVM=ON -DWITH_MPC=OFF
      -DWITH_OPENMP=ON -DWITH_SYMENGINE_THREAD_SAFE=ON)
  execute_process(COMMAND ${CMAKE_COMMAND} --build ${SYMENGINE_BUILD} --target
                          install)
endif()

find_package(SymEngine CONFIG REQUIRED PATHS ${SYMENGINE_STAGE})

find_package(PkgConfig REQUIRED)
pkg_check_modules(GMP REQUIRED IMPORTED_TARGET gmp)
pkg_check_modules(MPFR REQUIRED IMPORTED_TARGET mpfr)

FetchContent_Declare(
  cxxopts
  GIT_REPOSITORY https://github.com/jarro2783/cxxopts.git
  GIT_TAG v3.2.0)
FetchContent_MakeAvailable(cxxopts)
