
include_directories (${PROJECT_SOURCE_DIR}/pink)

add_library(pink
             lattribheight.c   
             llabelextrema.c   
             lwshedtopo.c      
             mccodimage.c      
             mccomptree.c      
             mcfahsalembier.c  
             mcimage.c         
             mcindic.c         
             mclifo.c          
             mcunionfind.c)

set (CMAKE_C_FLAGS "-fpic")
target_link_libraries(pink m)
