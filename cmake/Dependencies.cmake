include(FetchContent)

function(mbd_setup_dependencies)

  find_package(GiNaC REQUIRED)

  FetchContent_Declare(
    cxxopts
    GIT_REPOSITORY https://github.com/jarro2783/cxxopts.git
    GIT_TAG v3.2.0)
  FetchContent_MakeAvailable(cxxopts)

endfunction()
