cmake_minimum_required(VERSION 3.12 FATAL_ERROR)

set(MASTER_PROJECT OFF)
if(NOT DEFINED PROJECT_NAME)
  set(MASTER_PROJECT ON)
endif()

option(OPTIONS_WITH_EXAMPLES "Enable examples." ${MASTER_PROJECT})
option(OPTIONS_WITH_INSTALL "Enable installation." ${MASTER_PROJECT})
option(BUILD_SHARED_LIBS "Build shared libraries." OFF)

project(options VERSION 1.0 LANGUAGES CXX)

if(NOT CMAKE_BUILD_TYPE)
    message(STATUS "Setting build type to Release as none was specified.")
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY VALUE "Release")
endif()

if(NOT TARGET options)
    add_library(options ${options_SOURCE_DIR}/src/options.cpp)
    target_include_directories(options PUBLIC 
      $<BUILD_INTERFACE:${options_SOURCE_DIR}/src>
      )
    target_compile_features(options PUBLIC cxx_std_11)
endif()

if(OPTIONS_WITH_INSTALL)
    include(GNUInstallDirs)
    include(CMakePackageConfigHelpers)
    
    install(TARGETS options
            EXPORT options_targets 
            LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
            ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
            INCLUDES DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}")
    
    install(EXPORT options_targets
            FILE optionsConfig.cmake
            NAMESPACE options::
            DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/options")
    
    install(FILES src/options.hpp
            DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}")
    
    write_basic_package_version_file(
      "${options_BINARY_DIR}/optionsConfigVersion.cmake"
      VERSION ${options_VERSION}
      COMPATIBILITY SameMajorVersion)
    
    install(FILES "${options_BINARY_DIR}/optionsConfigVersion.cmake"
            DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/options")
endif()

if(OPTIONS_WITH_EXAMPLES)
    add_subdirectory(examples)
endif()
