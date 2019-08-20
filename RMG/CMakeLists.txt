list(APPEND RMGLIBS RMG_Common)
list(APPEND RMGLIBS RMG_Spin)
list(APPEND RMGLIBS RMG_Subdiag)
list(APPEND RMGLIBS zfp)
list(APPEND RMGLIBS symspg_static)
list(APPEND RMGLIBS RmgTddft)
list(APPEND RMGLIBS IO_Tddft)
list(APPEND RMGLIBS eldyn_mod)
list(APPEND RMGLIBS Eldyn)
if(RMG_GPU_ENABLED)
    list(APPEND RMGLIBS Gpufuncs)
endif()
if(QMCPACK_SUPPORT)
    list(APPEND RMGLIBS Interfaces)
endif()
if(QMCPACK_SUPPORT)
    list(APPEND RMGLIBS Interfaces)
endif()


add_executable(${rmg} Common/Main.cpp ${PROJECT_SOURCE_DIR}/Headers/main.h)
set_target_properties(${rmg} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")

target_include_directories(${rmg} PRIVATE "${PROJECT_SOURCE_DIR}/RMG/Headers/")


# There must be a better way to resolve circular references than this
foreach(lib ${RMGLIBS})
    target_link_libraries(${rmg} ${lib})
endforeach(lib ${RMGLIBS})
foreach(lib ${RMGLIBS})
    target_link_libraries(${rmg} ${lib})
endforeach(lib ${RMGLIBS})
foreach(lib ${RMGLIBS})
    target_link_libraries(${rmg} ${lib})
endforeach(lib ${RMGLIBS})
foreach(lib ${RMGLIBS})
    target_link_libraries(${rmg} ${lib})
endforeach(lib ${RMGLIBS})

install(TARGETS ${rmg} DESTINATION $PROJECT_BINARY_DIR)

