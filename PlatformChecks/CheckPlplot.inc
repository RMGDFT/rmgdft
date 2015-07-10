set(CMAKE_REQUIRED_INCLUDES "${PLplot_INCLUDE_DIR};${PLplot_INCLUDE_DIR}/plplot")
set(CMAKE_REQUIRED_LIBRARIES "${PLplot_LIBRARY}")
CHECK_C_SOURCE_COMPILES("#include <plplot.h>
                       int main(int argc, char **argv) {
                       plinit();
                     }" USE_PLPLOT_LIBS_INTERNAL)
