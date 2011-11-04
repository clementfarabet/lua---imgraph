
include_directories (${PROJECT_SOURCE_DIR}/PWcut ${PROJECT_SOURCE_DIR}/pink ${PROJECT_SOURCE_DIR}/graphcuts)

add_library(PWcut
	union_find.cpp 
	csparse.c
	MSF_RW.cpp 
	random_walker.cpp
	powerwatershed.cpp)

set (CMAKE_CXX_FLAGS "-fpic -g -DUNIXIO")
target_link_libraries(PWcut)
