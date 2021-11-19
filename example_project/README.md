# cpp-kde-fscr Example Project

This is an example on how to use **cpp-kde-fscr** library.

Make sure that cpp-kde-fscr library is already installed.

## Usage
``` bash

cmake -S . -B build -DCMAKE_PREFIX_PATH='/cpp-kde-fscr/installation/location'
cmake --build build
```

In `CMakeLists.txt`, we just need to add
``` cmake
find_package(cpp-kde-fscr CONFIG REQUIRED)
if (cpp-kde-fscr_FOUND)
  message(STATUS "cpp-kde-fscr library is found")
endif()

target_link_libraries(pdf PRIVATE cpp-kde-fscr::cpp-kde-fscr) # installed cpp-kde-fscr include/ path is automatically added
``` 
