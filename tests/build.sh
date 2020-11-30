#!/bin/bash
g++ -std=c++11 -I../include -o test_RectGrid test_RectGrid.cpp ../src/rectGrid.cpp ../src/fitsInterface.cpp -lCCfits -lcfitsio
