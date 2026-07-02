cmake_minimum_required(VERSION 3.16)

project(preshaires VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_library(preshaires-core
  src/core/pair_production.cpp
)
add_library(preshaires::core ALIAS preshaires-core)
target_include_directories(preshaires-core
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)
target_compile_options(preshaires-core PRIVATE -Wall -Wextra -Wpedantic)

add_library(preshaires-aires-adapter
  adapters/aires/aires_adapter.cpp
)
target_link_libraries(preshaires-aires-adapter PRIVATE preshaires::core)
target_compile_options(preshaires-aires-adapter PRIVATE -Wall -Wextra -Wpedantic)

add_executable(preshaires-cli
  apps/preshaires/main.cpp
)
set_target_properties(preshaires-cli PROPERTIES OUTPUT_NAME preshaires)
target_link_libraries(preshaires-cli PRIVATE preshaires::core)
target_compile_options(preshaires-cli PRIVATE -Wall -Wextra -Wpedantic)

enable_testing()

add_executable(preshaires-tests
  tests/unit/test_pair_production.cpp
)
target_link_libraries(preshaires-tests PRIVATE preshaires::core)
target_compile_options(preshaires-tests PRIVATE -Wall -Wextra -Wpedantic)

add_test(NAME preshaires-tests COMMAND preshaires-tests)
