cmake_minimum_required(VERSION 3.20)
project(benchmarks)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the type of build: Debug, Release, RelWithDebInfo, MinSizeRel"
      FORCE)
endif()

set(CMAKE_EXPORT_COMPILE_COMMANDS
    ON CACHE BOOL "Generate compile_commands.json")

# Add the ORC directory to CMAKE_PREFIX_PATH to help find the ORC package configuration
# set(CMAKE_PREFIX_PATH "/usr/lib/cmake/orc" ${CMAKE_PREFIX_PATH})
# set(ORC_DIR "/usr/lib/cmake/orc")

# Find the ORC package
find_package(ORC REQUIRED )

# Find ROOT package
find_package(ROOT REQUIRED
  COMPONENTS
    RIO
    Net
    ROOTNTuple
    Hist
    Graf
    Gpad
)

# Add the executable and link the libraries
add_executable(benchmarks benchmarks.cpp)

target_link_libraries(benchmarks
  PRIVATE
    orc::orc
    ROOT::ROOTNTuple
    ROOT::Hist
    ROOT::Graf
    ROOT::Gpad
)
