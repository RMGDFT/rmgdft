set(RMGDFT_PATH /home/luw/rmgdft)
include(ProcessorCount)
ProcessorCount(N)
if(NOT N EQUAL 0)
  set(CTEST_BUILD_FLAGS -j${N})
endif()

FILE(COPY ${RMGDFT_PATH}/CTestConfig.cmake  DESTINATION ${RMGDFT_PATH}/build)
include(${RMGDFT_PATH}/CTestConfig.cmake)
set(CTEST_CUSTOM_MAXIMUM_NUMBER_OF_ERRORS           "200" )
set(CTEST_CUSTOM_MAXIMUM_NUMBER_OF_WARNINGS         "500" )
set(CTEST_CUSTOM_MAXIMUM_PASSED_TEST_OUTPUT_SIZE    "104857600") # 100 MB
set(CTEST_CUSTOM_COVERAGE_EXCLUDE                   "")
# either customize these directly or write as CMake Template
# and use configure_file(... @ONLY) with CMake
set(CTEST_SOURCE_DIRECTORY   "${RMGDFT_PATH}/build")
set(CTEST_BINARY_DIRECTORY   "${RMGDFT_PATH}/build")
# build options
set(OPTION_BUILD             "-j8")
# define generator (optional), e.g. default to 'Unix Makefiles' on UNIX, Visual Studio on Windows
set(CTEST_GENERATOR          "Unix Makefiles")
# submit under Continuous, Nightly (default), Experimental
#set(CTEST_MODEL              "Nightly")
set(CTEST_MODEL              "Continuous")
# define how to checkout code, e.g. copy a directory, git pull, svn co, etc.
set(CTEST_CHECKOUT_COMMAND   "git pull")
# define how to update (optional), e.g. git checkout <git-branch>
#set(CTEST_UPDATE_COMMAND     "...")
# define how to configure (e.g. cmake -DCMAKE_INSTALL_PREFIX=...)
set(CTEST_CONFIGURE_COMMAND  "cmake ..")
# the name of the build
set(CTEST_BUILD_COMMAND      "make -j8 rmg-on-cpu rmg-cpu")
set(CTEST_BUILD_NAME         "rmg")
# how to build
# default max time each tests can run (in seconds)
set(CTEST_TIMEOUT            "7200")
# locale to English
set(ENV{LC_MESSAGES}         "en_EN")
ctest_start             (${CTEST_MODEL} TRACK ${CTEST_MODEL})
ctest_configure         (BUILD ${CTEST_BINARY_DIRECTORY} RETURN_VALUE ret_con)
ctest_build             (BUILD ${CTEST_BINARY_DIRECTORY} RETURN_VALUE ret_bld)

    # add as desired
ctest_test              (BUILD ${CTEST_BINARY_DIRECTORY} RETURN_VALUE ret_tst)
#    ctest_memcheck          (BUILD ${CTEST_BINARY_DIRECTORY} RETURN_VALUE ret_mem)
#    ctest_coverage          (BUILD ${CTEST_BINARY_DIRECTORY} RETURN_VALUE ret_cov)

    # attach build notes if desired, e.g. performance info, output files from tests
# list(APPEND CTEST_NOTES_FILES "/file/to/attach/as/build-note")

# standard submit
ctest_submit(RETURN_VALUE ret_sub)
# if dashboard requires a token that restricts who can submit to dashboard
#ctest_submit(RETURN_VALUE ret_sub HTTPHEADER "Authorization: Bearer ${CTEST_TOKEN}")
