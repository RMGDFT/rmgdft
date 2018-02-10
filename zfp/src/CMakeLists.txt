#if(UNIX)
#  target_link_libraries(zfp PUBLIC m)
#endif()
#
#if(WIN32)
#  # Define ZFP_SOURCE when compiling libzfp to export symbols to Windows DLL
#  list(APPEND zfp_defs ZFP_SOURCE)
#endif()
#
#target_compile_definitions(zfp PRIVATE ${zfp_defs})
#
#target_include_directories(zfp
#  PUBLIC    ${ZFP_SOURCE_DIR}/include
#  INTERFACE ${ZFP_SOURCE_DIR}/array)
#
#set_property(TARGET zfp PROPERTY VERSION ${ZFP_VERSION})
#set_property(TARGET zfp PROPERTY SOVERSION 0)

