cmake_minimum_required(VERSION 3.14)

include(FetchContent) # https://cmake.org/cmake/help/v3.14/module/FetchContent.html
FetchContent_Declare(
  googletest
  # Specify the commit you depend on and update it regularly.
  URL https://github.com/google/googletest/archive/609281088cfefc76f9d0ce82e1ff6c30cc3591e5.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

include_directories(${PROJECT_SOURCE_DIR}/include)

# List all files containing tests
set(TESTFILES
  kde-test.cpp
)
set(TEST_MAIN unit_tests)

# --------------------------------------------------------------------------------
#                                   Make Tests
# --------------------------------------------------------------------------------
add_executable(${TEST_MAIN} ${TESTFILES})
target_link_libraries(
  ${TEST_MAIN}
  PRIVATE gtest_main ${PROJECT_NAME}
)

add_custom_target(tests)
add_dependencies(tests ${TEST_MAIN})

set_target_properties(${TEST_MAIN} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})
set_target_properties(
  ${TEST_MAIN}
    PROPERTIES
      CXX_STANDARD 11
      CXX_STANDARD_REQUIRED YES
)

include(GoogleTest) # https://cmake.org/cmake/help/v3.14/module/GoogleTest.html
gtest_discover_tests(${TEST_MAIN})

include(CTest) # https://cmake.org/cmake/help/v3.14/module/CTest.html
