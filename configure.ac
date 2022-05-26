AC_INIT([libvkl], [0.1], [gvernard@astro.rug.nl])
AC_CONFIG_AUX_DIR(aux-dist)
AM_INIT_AUTOMAKE([foreign subdir-objects])
AC_CONFIG_MACRO_DIRS([m4])

AC_PROG_CXX
AC_PROG_F77

LT_INIT([disable-static])
AC_FUNC_MALLOC
AC_CHECK_FUNCS([floor pow sqrt])



# adding the lib to the files to link
##########################################################################
# adding the XXX library (e.g. with static name 'libXXX.a')
##########################################################################
LIBS="-lgfortran -lCCfits -lcfitsio -lgmp -lCGAL $LIBS"

# introduce the optional configure parameter for a non-standard install prefix of XXX
AC_ARG_WITH([external],
    [AS_HELP_STRING([--with-external=prefix],[try this for a non-standard install prefix of all the external libraries])],
    [EXT_PATHSET=1],
    [EXT_PATHSET=0])
	
# if optional parameter used, extend path flags for compliler and linker
if test $EXT_PATHSET = 1 ; then
    # extend the compiler and linker flags according to the path set
    AM_CXXFLAGS="-I${with_external}cfitsio/include -I${with_external}CCfits/include $AM_CXXFLAGS"
    AM_LDFLAGS="$AM_LDFLAGS -L${with_external}cfitsio/lib -L${with_external}CCfits/lib"
fi

# distribute the changed variables among the Makefiles
AC_SUBST([LIBS])
AC_SUBST([AM_CXXFLAGS])
AC_SUBST([AM_LDFLAGS])



AC_CONFIG_FILES([Makefile])

AC_OUTPUT
