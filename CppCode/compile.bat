c++ -O3 -g -Wall -shared -std=c++11 -fPIC $(python3-config --includes) -I. -I../extern/pybind11/include ht.cpp htbind.cpp -o ht$(python3-config --extension-suffix)
