include(ExternalProject)

function(mbd_setup_dependencies)
  set(CUSTOM_CFLAGS "-O3 -march=native -mtune=native")
  set(CUSTOM_CXXFLAGS "-O3 -march=native -mtune=native")
  set(CUSTOM_LDFLAGS "-flto")

  set(CLN_SOURCE_DIR ${CMAKE_SOURCE_DIR}/libs/cln)

  add_subdirectory(libs/ginac)

endfunction()
