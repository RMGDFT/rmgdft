# options

This is a small library for command-line options. It does not require any additional dependencies.

## Building and Installing

Assuming that you want to use the `gcc 8` compiler and `OpenMP`, you can build the project as follows:
```bash
# clone the repo
git clone https://github.com/kabicm/options
cd options
mkdir build
cd build

# build
CC=gcc-8 CXX=g++-8 cmake ..

# compile
make -j 4
```

## Example

The library is very simple and easy to use:
```cpp
// first initialize the library with command line arguments
options::initialize(argc, argv);

// default value will be used if this option is not specified
int default_value = 10;
int m = options::next_int("-m", "--m_value", "Description of the value.", default_value);
```
There is also a small example that reads 2 ints from the command-line and outputs them. The example can be run within build folder as:
```bash
./examples/simple -m 10 -n 10
```

## Author
Marko Kabic (marko.kabic@cscs.ch)
