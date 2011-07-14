
package = "imgraph"
version = "1.0-1"

source = {
   url = "imgraph-1.0-1.tgz"
}

description = {
   summary = "A package to deal with graphs on images.",
   detailed = [[
            This package provides standard functions to
            create and manipulate edge-weighted graphs 
            of images.
   ]],
   homepage = "",
   license = "MIT/X11" -- or whatever you like
}

dependencies = {
   "lua >= 5.1",
   "torch",
   "sys",
   "xlua",
   "image"
}

build = {
   type = "cmake",

   cmake = [[
         cmake_minimum_required(VERSION 2.8)

         set (CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR})

         # infer path for Torch7
         string (REGEX REPLACE "(.*)lib/luarocks/rocks.*" "\\1" TORCH_PREFIX "${CMAKE_INSTALL_PREFIX}" )
         message (STATUS "Found Torch7, installed in: " ${TORCH_PREFIX})

         find_package (Torch REQUIRED)

         set (CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

         include_directories (${TORCH_INCLUDE_DIR} ${PROJECT_SOURCE_DIR})
         add_library (imgraph SHARED init.c)
         link_directories (${TORCH_LIBRARY_DIR})
         target_link_libraries (imgraph ${TORCH_LIBRARIES})

         install_files(/lua/imgraph init.lua)
         install_targets(/lib imgraph)
   ]],

   variables = {
      CMAKE_INSTALL_PREFIX = "$(PREFIX)"
   }
}
