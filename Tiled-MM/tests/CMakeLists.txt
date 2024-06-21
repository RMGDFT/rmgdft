################
#  Build test  #
################

add_executable(test-multiply "test-multiply.cpp")
target_link_libraries(test-multiply PRIVATE Tiled-MM cxxopts::cxxopts)

################
#  unit gtests #
################
add_test(NAME square-small COMMAND test-multiply -m 1000 -n 1000 -k 1000)
add_test(NAME square-large COMMAND test-multiply -m 10000 -n 10000 -k 10000)
add_test(NAME non-square-small COMMAND test-multiply -m 1234 -n 4567 -k 1357)
add_test(NAME non-square-large COMMAND test-multiply -m 12345 -n 23456 -k 67891)
