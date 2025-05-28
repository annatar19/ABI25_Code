cmake_minimum_required(VERSION 3.20)
project(lhcb_rntuple_orc)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the type of build: Debug, Release, RelWithDebInfo, MinSizeRel"
      FORCE)
endif()

set(CMAKE_EXPORT_COMPILE_COMMANDS
    ON CACHE BOOL "Generate compile_commands.json")

set(ORC_DIR "/home/annatar/repos/orc/build/lib/cmake/orc")
list(APPEND CMAKE_PREFIX_PATH "${ORC_DIR}")

find_package(ORC CONFIG REQUIRED
    COMPONENTS orc
)


set(ROOT_DIR "/home/annatar/studie/abi25/root/root-build"
    CACHE PATH "Path to the ROOT installation")

find_package(ROOT CONFIG REQUIRED
  COMPONENTS
    RIO
    Net
    ROOTNTuple
    Hist      # for TH1D, TH2D, …
    Graf      # for TCanvas
    Gpad      # for TPad
)

add_executable(lhcb_rntuple_orc lhcb_rntuple_orc.cpp)

target_link_libraries(lhcb_rntuple_orc
  PRIVATE
    orc::orc
    ROOT::ROOTNTuple
    ROOT::Hist
    ROOT::Graf
    ROOT::Gpad
)
