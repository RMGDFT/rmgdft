list(APPEND RMGLIBS NEGF_Common)
list(APPEND RMGLIBS NEGF_Input)
list(APPEND RMGLIBS NEGF_Pois)
list(APPEND RMGLIBS NEGF_Lcurrent)
list(APPEND RMGLIBS NEGF_Green)
list(APPEND RMGLIBS ON_NEGF_share)

add_executable(${rmg-negf} EXCLUDE_FROM_ALL Common/Main.cpp ${PROJECT_SOURCE_DIR}/Headers/main.h)
set_target_properties(${rmg} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")

target_include_directories(${rmg-negf} PRIVATE "${PROJECT_SOURCE_DIR}/NEGF/Headers/")


# There must be a better way to resolve circular references than this
foreach(lib ${RMGLIBS})
    target_link_libraries(${rmg-negf} ${lib})
endforeach(lib ${RMGLIBS})
foreach(lib ${RMGLIBS})
    target_link_libraries(${rmg-negf} ${lib})
endforeach(lib ${RMGLIBS})
foreach(lib ${RMGLIBS})
    target_link_libraries(${rmg-negf} ${lib})
endforeach(lib ${RMGLIBS})
foreach(lib ${RMGLIBS})
    target_link_libraries(${rmg-negf} ${lib})
endforeach(lib ${RMGLIBS})

#install(TARGETS ${rmg-negf} DESTINATION bin)

