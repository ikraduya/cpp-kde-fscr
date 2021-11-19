cmake_minimum_required(VERSION 3.14)

project(Example_cpp-kde-fscr)

find_package(cpp-kde-fscr CONFIG REQUIRED)
if (cpp-kde-fscr_FOUND)
  message(STATUS "cpp-kde-fscr library is found")
endif()

# An executable for the file(s)
add_executable(pdf pdf.cpp)
target_link_libraries(pdf PRIVATE cpp-kde-fscr::cpp-kde-fscr)

# Set the properties required, e.g. C++ standard used. Apply to library and examples
set_target_properties(
  pdf
  PROPERTIES
    CXX_STANDARD 11
    CXX_STANDARD_REQUIRED YES
)