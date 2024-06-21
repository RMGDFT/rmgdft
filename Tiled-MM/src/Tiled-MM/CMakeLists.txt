add_library(Tiled-MM gpu_context.cpp
                     mm_handle.cpp
                     tile_dim.cpp
                     tiled_matrix.cpp
                     tiled_mm.cpp
                     tile_coord.cpp)

target_link_libraries(Tiled-MM PUBLIC
  $<TARGET_NAME_IF_EXISTS:tiledmm::rocm>
  $<TARGET_NAME_IF_EXISTS:tiledmm::cuda>
)

message("TILEDMM_EXTERNAL_LIBRARIES ${TILEDMM_EXTERNAL_LIBRARIES}")

target_include_directories(Tiled-MM PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
    )

target_compile_features(Tiled-MM PUBLIC cxx_std_14)

target_compile_definitions(Tiled-MM PUBLIC ${TILEDMM_DEFINITIONS})

install(TARGETS Tiled-MM
  EXPORT Tiled-MM_targets
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  INCLUDES DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}")
