// Local
#include <options.hpp>

int main( int argc, char **argv ) {
    options::initialize(argc, argv);

    auto m = options::next_int("-m", "--rows", "Number of rows.", 10);
    auto n = options::next_int("-n", "--cols", "Number of cols.", 12);

    std::cout << "Number of rows: " << m << std::endl;
    std::cout << "Number of cols: " << n << std::endl;

    return 0;
}
