cmake_minimum_required (VERSION 3.5)
project(pyClustering)
find_package(pybind11 REQUIRED)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories(/Users/cdesio/anaconda3/envs/clustering/lib/python3.6/site-packages/pybind11/include/pybind11/)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)

set(pybind11_DIR /Users/cdesio/anaconda3/envs/clustering/share/cmake/pybind11/)

include(cmake/CheckGit.cmake)
CheckGitSetup()

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

file(GLOB MAC_FILES 
${PROJECT_SOURCE_DIR}/runSingle.py
${PROJECT_SOURCE_DIR}/sugarPosTest.csv
${PROJECT_SOURCE_DIR}/test.csv
${PROJECT_SOURCE_DIR}/testFile.root
)

foreach(_script ${MAC_FILES})
  configure_file(
    ${_script}
    ${PROJECT_BINARY_DIR}/.
    COPYONLY
    )
endforeach()

pybind11_add_module(clustering pyClustering.cc  ${sources} ${headers})
