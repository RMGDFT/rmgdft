add_library(zfp
  src/zfp.c
  src/bitstream.c
  src/traitsf.h src/traitsd.h src/block1.h src/block2.h src/block3.h
  src/encode1f.c src/encode1d.c src/encode1i.c src/encode1l.c
  src/decode1f.c src/decode1d.c src/decode1i.c src/decode1l.c
  src/encode2f.c src/encode2d.c src/encode2i.c src/encode2l.c
  src/decode2f.c src/decode2d.c src/decode2i.c src/decode2l.c
  src/encode3f.c src/encode3d.c src/encode3i.c src/encode3l.c
  src/decode3f.c src/decode3d.c src/decode3i.c src/decode3l.c
)


#------------------------------------------------------------------------------#
# Parse version number from zfp.h
#------------------------------------------------------------------------------#
file(READ ${CMAKE_CURRENT_SOURCE_DIR}/include/zfp.h _zfp_h_contents)
string(REGEX REPLACE ".*#define[ \t]+ZFP_VERSION_MAJOR[ \t]+([0-9]+).*"
     "\\1" ZFP_VERSION_MAJOR ${_zfp_h_contents})
string(REGEX REPLACE ".*#define[ \t]+ZFP_VERSION_MINOR[ \t]+([0-9]+).*"
    "\\1" ZFP_VERSION_MINOR ${_zfp_h_contents})
string(REGEX REPLACE ".*#define[ \t]+ZFP_VERSION_PATCH[ \t]+([0-9]+).*"
    "\\1" ZFP_VERSION_PATCH ${_zfp_h_contents})
set(ZFP_VERSION
  "${ZFP_VERSION_MAJOR}.${ZFP_VERSION_MINOR}.${ZFP_VERSION_PATCH}")

#project(ZFP VERSION ${ZFP_VERSION})

#------------------------------------------------------------------------------#
# Top level options
#------------------------------------------------------------------------------#

# Windows (Visual Studio) specific options
if(MSVC)
  # Use this to get a usable export library when building a DLL on Windows
  set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)

  # Silence extraneous Visual Studio specific warnings
  add_definitions(-D_CRT_SECURE_NO_WARNINGS -D_SCL_SECURE_NO_WARNINGS /wd4146 /wd4305)
endif()

# Suggest C99
if(NOT CMAKE_C_STANDARD)
  set(CMAKE_C_STANDARD 99)
endif()

if(MSVC OR MINGW)
  set(CMAKE_C_STANDARD 90)
endif()

message(STATUS "Compiling with C standard: ${CMAKE_C_STANDARD}")

# Suggest C++98
if(NOT CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 98)
endif()
message(STATUS "Compiling with C++ standard: ${CMAKE_CXX_STANDARD}")

include(CMakeDependentOption)


set(ZFP_BIT_STREAM_WORD_SIZE 64 CACHE STRING
  "Use smaller bit stream word type for finer rate granularity")
set_property(CACHE ZFP_BIT_STREAM_WORD_SIZE PROPERTY STRINGS "8;16;32;64")

option(ZFP_WITH_BIT_STREAM_STRIDED
  "Enable strided access for progressive zfp streams" OFF)

option(ZFP_WITH_ALIGNED_ALLOC "Enabled aligned memory allocation" OFF)

option(ZFP_WITH_CACHE_TWOWAY "Use two-way skew-associative cache" OFF)

option(ZFP_WITH_CACHE_FAST_HASH
  "Use a faster but more collision prone hash function" OFF)

option(ZFP_WITH_CACHE_PROFILE "Count cache misses" OFF)

# Handle compile-time macros

if((DEFINED ZFP_INT64) AND (DEFINED ZFP_INT64_SUFFIX))
  list(APPEND zfp_defs ZFP_INT64=${ZFP_INT64})
  list(APPEND zfp_defs ZFP_INT64_SUFFIX=${ZFP_INT64_SUFFIX})
endif()

if((DEFINED ZFP_UINT64) AND (DEFINED ZFP_UINT64_SUFFIX))
  list(APPEND zfp_defs ZFP_UINT64=${ZFP_UINT64})
  list(APPEND zfp_defs ZFP_UINT64_SUFFIX=${ZFP_UINT64_SUFFIX})
endif()

if(NOT (ZFP_BIT_STREAM_WORD_SIZE EQUAL 64))
  list(APPEND zfp_defs BIT_STREAM_WORD_TYPE=uint${ZFP_BIT_STREAM_WORD_SIZE})
endif()

if(ZFP_WITH_BIT_STREAM_STRIDED)
  list(APPEND zfp_defs BIT_STREAM_STRIDED)
endif()

if(ZFP_WITH_ALIGNED_ALLOC)
  list(APPEND zfp_defs ZFP_ALIGNED_ALLOC)
endif()

if(ZFP_WITH_CACHE_TWOWAY)
  list(APPEND zfp_defs ZFP_CACHE_TWOWAY)
endif()

if(ZFP_WITH_CACHE_FAST_HASH)
  list(APPEND zfp_defs ZFP_CACHE_FAST_HASH)
endif()

if(ZFP_WITH_CACHE_PROFILE)
  list(APPEND zfp_defs ZFP_CACHE_PROFILE)
endif()

#------------------------------------------------------------------------------#
# Add source code
#------------------------------------------------------------------------------#
#add_subdirectory(src)

#------------------------------------------------------------------------------#
# Header install
#------------------------------------------------------------------------------#

#install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
#install(DIRECTORY array/   DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

