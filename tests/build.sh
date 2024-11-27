#!/bin/bash
#g++ -std=c++11 -I../include -o test_RectGrid test_RectGrid.cpp ../src/rectGrid.cpp ../src/fitsInterface.cpp -lCCfits -lcfitsio

mylibs="/home/giorgos/myLibraries/libraries/"
g++ -std=c++11 -o test_multipole test_multipole.cpp -I${mylibs}vkl_lib/include -I${mylibs}jsoncpp/include -L${mylibs}/vkl_lib/lib -L${mylibs}jsoncpp/lib -Wl,-rpath=${mylibs}vkl_lib/lib -Wl,-rpath=${mylibs}jsoncpp/lib -ljsoncpp -lvkl
