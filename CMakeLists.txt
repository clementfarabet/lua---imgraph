
CMAKE_MINIMUM_REQUIRED(VERSION 2.6 FATAL_ERROR)
CMAKE_POLICY(VERSION 2.6)
IF(LUAROCKS_PREFIX)
    MESSAGE(STATUS "Installing Torch through Luarocks")
    STRING(REGEX REPLACE "(.*)lib/luarocks/rocks.*" "\\1" CMAKE_INSTALL_PREFIX  "${LUAROCKS_PREFIX}")
    MESSAGE(STATUS "Prefix inferred from Luarocks: ${CMAKE_INSTALL_PREFIX}")
ENDIF()
FIND_PACKAGE(Torch REQUIRED)
FIND_PACKAGE(Boost REQUIRED)

add_subdirectory (pink)
add_subdirectory (graphcuts)
add_subdirectory (PWcut)

include_directories (${PROJECT_SOURCE_DIR} ${PROJECT_SOURCE_DIR}/pink
${PROJECT_SOURCE_DIR}/graphcuts ${PROJECT_SOURCE_DIR}/PWcut ${Boost_INCLUDE_DIR})

SET(src init.cpp)
SET(luasrc init.lua MalisCriterion.lua)

ADD_TORCH_PACKAGE(imgraph "${src}" "${luasrc}" "Image Processing")
SET (CMAKE_CXX_FLAGS "-fpic")
TARGET_LINK_LIBRARIES(imgraph luaT TH graphcuts pink PWcut)
