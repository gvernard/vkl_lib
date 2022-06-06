AC_INIT([libvkl], [0.1], [gvernard@astro.rug.nl])
AC_LANG(C++)
AC_CONFIG_AUX_DIR(aux-dist)
AM_INIT_AUTOMAKE([foreign subdir-objects])
AC_CONFIG_MACRO_DIRS([m4])
AC_PREFIX_DEFAULT($PWD)

AC_PROG_CXX
AC_PROG_F77

LT_INIT([disable-static])


ac_new_CPPFLAGS="-Iinclude"
ac_new_LDFLAGS=""
ac_new_LIBS=""

AC_ARG_WITH([cfitsio],
    [AS_HELP_STRING([--with-cfitsio=prefix],[give the path to the cfitsio library and header files.])],
    [
	ac_new_CPPFLAGS+=" -I${withval}/include"
	ac_new_LDFLAGS+=" -L${withval}/lib -Wl,-rpath -Wl,${withval}/lib"
	ac_new_LIBS+=" -lcfitsio"
    ],
    [
	AC_CHECK_LIB(cfitsio,main,[],[AC_MSG_ERROR("you need to have the cfitsio library installed.")])
	AC_CHECK_HEADERS(fitsio.h,[],[AC_MSG_ERROR("you need to have the cfitsio header files installed.")])
    ])

AC_ARG_WITH([CCfits],
    [AS_HELP_STRING([--with-CCfits=prefix],[give the path to the CCfits library and header files.])],
    [
	ac_new_CPPFLAGS+=" -I${withval}/include"
	ac_new_LDFLAGS+=" -L${withval}/lib -Wl,-rpath -Wl,${withval}/lib"
	ac_new_LIBS+=" -lCCfits"
    ],
    [
	AC_CHECK_LIB(CCfits,main,[],[AC_MSG_ERROR("you need to have the CCfits library installed.")],[-lcfitsio])
	AC_CHECK_HEADERS(CCfits/CCfits.h,[],[AC_MSG_ERROR("you need to have the CCfits header files installed.")])
    ])

AC_ARG_WITH([gmp],
    [AS_HELP_STRING([--with-gmp=prefix],[give the path to the gmp library and header files.])],
    [
	ac_new_CPPFLAGS+=" -I${withval}/include"
	ac_new_LDFLAGS+=" -L${withval}/lib -Wl,-rpath -Wl,${withval}/lib"
	ac_new_LIBS+=" -lgmp"
    ],
    [
	AC_CHECK_LIB(gmp,main,[],[AC_MSG_ERROR("you need to have the gmp library installed.")])
	AC_CHECK_HEADERS(gmp.h,[],[AC_MSG_ERROR("you need to have the gmp header files installed.")])
    ])

AC_ARG_WITH([CGAL],
    [AS_HELP_STRING([--with-CGAL=prefix],[give the path to the CGAL library and header files.])],
    [
	ac_new_CPPFLAGS+=" -I${withval}/include"
	ac_new_LDFLAGS+=" -L${withval}/lib64 -Wl,-rpath -Wl,${withval}/lib"
	ac_new_LIBS+=" -lCGAL"
    ],
    [
	AC_CHECK_LIB(CGAL,main,[],[AC_MSG_ERROR("you need to have the CGAL library installed.")])
	#AC_CHECK_HEADERS(CGAL/CGAL.h,[],[AC_MSG_ERROR("you need to have the CGAL header files installed.")])
    ])

AC_ARG_WITH([jsoncpp],
    [AS_HELP_STRING([--with-jsoncpp=prefix],[give the path to the jsoncpp library and header files.])],
    [
	ac_new_CPPFLAGS+=" -I${withval}/include"
	ac_new_LDFLAGS+=" -L${withval}/lib64 -Wl,-rpath -Wl,${withval}/lib"
	ac_new_LIBS+=" -ljsoncpp"
    ],
    [
	AC_CHECK_LIB(jsoncpp,main,[],[AC_MSG_ERROR("you need to have the jsoncpp library installed.")])
	AC_CHECK_HEADERS(json/json.h,[],[AC_MSG_ERROR("you need to have the jsoncpp header files installed.")])
    ])

CPPFLAGS+="$ac_new_CPPFLAGS"
LDFLAGS+="$ac_new_LDFLAGS"
LIBS+="$ac_new_LIBS"


AC_SUBST([CPPFLAGS])
AC_SUBST([LDFLAGS])
AC_SUBST([LIBS])

AC_CONFIG_FILES([Makefile])

AC_OUTPUT
