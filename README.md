# VKLlib

This is a C++ library of functions to perform various lensing tasks, e.g. define a mass model, define a source, produce lensed images, etc.



## Getting Started
### Prerequisites

There is a number of third-party libraries required to install vkllib, namely: **cfitsio**, **CCfits**, **gmp**, **CGAL**, and **jsoncpp*.
In additon, **cmake** and **autotools** (the autoreconf command) are required to perform the installation. 


### Installing

To install vkllib one needs to call `autoreconf -i' and then the usual './configure, make, make install'.
If the required third-party libraries listed above are not installed in a standard system location, they will need to be explicitly provided.
This is possible via '--with-<library_name>=/path/to/library' options passed to the configure script.
An example call to the configure script could look like the following:

```
./configure --prefix=/path/to/installation/of/vkl_lib --with-cfitsio=/path/to/libraries/cfitsio --with-CCfits=/path/to/libraries/CCfits --with-gmp=/path/to/libraries/gmp --with-CGAL=/path/to/libraries/CGAL --with-jsoncpp=/path/to/libraries/jsoncpp
```


### Finalizing

For convenience, in order to be able to compile programs that use gerlumphpp without the need to explicitly specify the path to the headers and the libarry (e.g. the -I, -L, and -Wl,-rpath options to the g++ compiler and linker), one can define the following environment variables:

```
CPATH=$CPATH:/path/to/vkllib/include
LIBRARY_PATH=$LIBRARY_PATH:/path/to/installation/of/vkllib/lib
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/path/to/installation/of/vkllib/lib
```

For the bash shell, one can export these variables from within the .bashrc file.




## Authors

**Georgios Vernardos** ( [github](https://github.com/gvernard)  - [homepage](http://astronomy.swin.edu.au/~gvernard/) )



