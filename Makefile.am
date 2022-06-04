ACLOCAL_AMFLAGS=-I m4

lib_LTLIBRARIES=libvkl.la

libvkl_la_SOURCES= \
	src/baseSourcePlane.cpp \
	src/covKernels.cpp \
	src/fixedSource.cpp \
	src/lightProfile.cpp \
	src/rectGrid.cpp \
	src/constants.cpp \
	src/fitsInterface.cpp \
	src/imagePlane.cpp \
	src/jsonParsers.cpp \
	src/massModels.cpp \
	src/fastell.f

libvkl_la_CXXFLAGS = -std=c++11	-fPIC -g -frounding-math

libvkl_la_LIBADD = -lgfortran -lCCfits -lcfitsio -lgmp -lCGAL -ljsoncpp
