cmake_minimum_required(VERSION 3.20)
project(rntuple_scaling)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the type of build: Debug, Release, RelWithDebInfo, MinSizeRel"
      FORCE)
endif()

set(CMAKE_EXPORT_COMPILE_COMMANDS
    ON CACHE BOOL "Generate compile_commands.json")

set(ROOT_DIR "/home/annatar/studie/abi25/root/root-build"
    CACHE PATH "Path to the ROOT installation")

find_package(ROOT CONFIG REQUIRED
    COMPONENTS
      ROOTNTuple
)

file(GLOB SOURCES
     "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp"
)

add_executable(rntuple_scaling ${SOURCES})

target_link_libraries(rntuple_scaling
    PRIVATE
      ROOT::ROOTNTuple
)
