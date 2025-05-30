cmake_minimum_required(VERSION 3.10) 

project(HENS)

add_compile_options(-Wall -Wextra -Wpedantic)

set(HENS_MATH_CORE_SRC src/HENS_math_core.cpp)
set(HENS_MATRIX_SRC src/HENS_matrix.cpp)
set(MAIN_SRC src/HENS.cpp)

add_library(hens_math_core ${HENS_MATH_CORE_SRC})
add_library(hens_matrix ${HENS_MATRIX_SRC})

target_include_directories(hens_math_core PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_include_directories(hens_matrix PUBLIC ${PROJECT_SOURCE_DIR}/include)

add_executable(HENS ${MAIN_SRC})

target_link_libraries(hens_math_core hens_matrix)
target_link_libraries(HENS hens_math_core hens_matrix)

target_include_directories(HENS PUBLIC ${PROJECT_SOURCE_DIR}/include)