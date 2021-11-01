cmake_minimum_required(VERSION 3.13)
project(SVcompare)

set(CMAKE_CXX_STANDARD 11)

include_directories(htslib-1.13/include)

add_executable(compare-ins lib/ssw.c lib/ssw_cpp.cpp compare-ins.cpp)
add_executable(compare-del lib/ssw.c lib/ssw_cpp.cpp compare-del.cpp)

find_library(HTS_LIB hts PATHS htslib-1.13/lib NO_DEFAULT_PATH)
target_link_libraries(compare-del "${HTS_LIB}")
target_link_libraries(compare-ins "${HTS_LIB}")
