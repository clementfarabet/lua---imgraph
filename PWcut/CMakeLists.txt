
include_directories (${PROJECT_SOURCE_DIR}/PWcut ${PROJECT_SOURCE_DIR}/pink ${PROJECT_SOURCE_DIR}/graphcuts)

add_library(PWcut
	union_find.cpp 
	csparse.c
	MSF_RW.cpp 
	random_walker.cpp
	powerwatershed.cpp
	mclistechainee.c)

set (CMAKE_CXX_FLAGS "-fpic -g -DUNIXIO")
set (CMAKE_LINKER_FLAGS "-fpic")
target_link_libraries(PWcut)
