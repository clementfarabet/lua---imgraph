
include_directories (${PROJECT_SOURCE_DIR}/graphcuts)

add_library(graphcuts
             maxflow.cpp)

set (CMAKE_CXX_FLAGS "-fpic -g -DUNIXIO")
target_link_libraries(graphcuts)
