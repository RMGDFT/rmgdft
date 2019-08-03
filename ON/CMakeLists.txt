list(APPEND RMGLIBS RmgLib)
list(APPEND RMGLIBS ON_Common)
list(APPEND RMGLIBS ON_Input)
list(APPEND RMGLIBS ON_NEGF_share)
list(APPEND RMGLIBS MG)
list(APPEND RMGLIBS zfp)
list(APPEND RMGLIBS RMG_Spin)
list(APPEND RMGLIBS symspg_static)
list(APPEND RMGLIBS Force)
list(APPEND RMGLIBS Input)
list(APPEND RMGLIBS Misc)
list(APPEND RMGLIBS XC)
list(APPEND RMGLIBS Finite_diff)
list(APPEND RMGLIBS US_PP)
if(QMCPACK_SUPPORT)
    list(APPEND RMGLIBS Interfaces)
endif()
if(NOT CRAY_XK)
    list(APPEND RMGLIBS scalapack)
endif()

add_executable(${rmg-on} EXCLUDE_FROM_ALL Common/Main.cpp ${PROJECT_SOURCE_DIR}/Headers/main.h)
set_target_properties(${rmg-on} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")

target_include_directories(${rmg-on} PRIVATE "${PROJECT_SOURCE_DIR}/ON/Headers/")

if(NOT CRAY_XK)
    add_dependencies(${rmg-on} scalapack)
endif(NOT CRAY_XK)


# There must be a better way to resolve circular references than this
foreach(lib ${RMGLIBS})
    target_link_libraries(${rmg-on} ${lib})
endforeach(lib ${RMGLIBS})
foreach(lib ${RMGLIBS})
    target_link_libraries(${rmg-on} ${lib})
endforeach(lib ${RMGLIBS})
foreach(lib ${RMGLIBS})
    target_link_libraries(${rmg-on} ${lib})
endforeach(lib ${RMGLIBS})
foreach(lib ${RMGLIBS})
    target_link_libraries(${rmg-on} ${lib})
endforeach(lib ${RMGLIBS})


#install(TARGETS ${rmg-on} DESTINATION bin)

