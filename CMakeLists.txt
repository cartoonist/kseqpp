cmake_minimum_required(VERSION 3.10)
project(kseq++ VERSION 0.1.4 LANGUAGES CXX)

# Finding dependencies.
find_package(ZLIB REQUIRED)
find_package(BZip2 REQUIRED)
find_package(Threads REQUIRED)

# Force CMake to use `-I` instead of `-isystem`
set(CMAKE_NO_SYSTEM_FROM_IMPORTED TRUE)

# Defining the configure header file
configure_file(include/kseq++/config.hpp.in ${PROJECT_SOURCE_DIR}/include/kseq++/config.hpp @ONLY)
# Setting header files
file(GLOB HEADER_FILES "${PROJECT_SOURCE_DIR}/include/kseq++/*.h")
# Creating an INTERFACE library
add_library(kseq++ INTERFACE)
# Defining target sources
target_sources(kseq++ INTERFACE "$<BUILD_INTERFACE:${HEADER_FILES}>")
# Defining include directories
target_include_directories(kseq++
  INTERFACE $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
  SYSTEM INTERFACE $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/include>
  INTERFACE ZLIB::ZLIB
  INTERFACE BZip2::BZip2
  INTERFACE Threads::Threads)
# Defining link libraries
target_link_libraries(kseq++
  INTERFACE ZLIB::ZLIB
  INTERFACE BZip2::BZip2
  INTERFACE Threads::Threads)
# Use C++17
target_compile_features(kseq++ INTERFACE cxx_std_11)
# Defining `kseq++::kseq++` alias
add_library(kseq++::kseq++ ALIAS kseq++)
# Generating pkgconfig file
configure_file(kseq++.pc.in ${PROJECT_SOURCE_DIR}/kseq++.pc @ONLY)
# Install targets
install(DIRECTORY include/kseq++ DESTINATION include PATTERN "*.hpp.in" EXCLUDE)
install(FILES kseq++.pc DESTINATION lib/pkgconfig)
