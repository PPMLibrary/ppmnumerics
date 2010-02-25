#-------------------------------------------------------------------------
#  File         :  configure.ac
#-------------------------------------------------------------------------
#
#  Purpose      :  Compilation
#
#  Remarks      :
#
#  References   :
#
#  Revisions    : This autoconf was written by Justin Park  8.17.2009
#-------------------------------------------------------------------------
#  Parallel Particle Mesh Library (PPM)
#  Institute of Computational Science
#  ETH Zentrum, Hirschengraben 84
#  CH-8092 Zurich, Switzerland
#-------------------------------------------------------------------------

# Initialization
AC_PREREQ(2.61)
m4_pattern_allow([AS_TR_SH])

AC_INIT(PPM, 1.0.2)
#AC_CONFIG_SRCDIR([src/ppm.h])

# Check architecture
AC_CHECK_SIZEOF([long double])

# Checks for basic compiler and flags
AC_PROG_CXX
AC_PROG_CPP
AC_PROG_MAKE_SET
AC_PROG_RANLIB
AC_PROG_FC
AC_FC_FREEFORM
AC_FC_LIBRARY_LDFLAGS

# for now let's use the C compiler to check the vendor we're using here
# this of course means that we HAVE to set all compiler flags to compilers of
# the same vendor (which is anyway a good idea)
AX_COMPILER_VENDOR

# Remove existing preprocessing flags
touch ./src/ppm_define.h
rm ./src/ppm_define.h

# Checks for libraries and set preprocessing flags

# this test should be omitted if the compiler is intel, they provide their own
# math implementation
if test "x${ax_cv_c_compiler_vendor}" != xintel; then
  AC_CHECK_LIB([m],[sqrt])
fi

# Make sure that the libraries are compiled with Fortran
AC_LANG([Fortran])

AX_F90_MODULE_EXTENSION
if test x$ax_cv_f90_modext = xunknown ; then
    AC_MSG_ERROR([unable to find f90 modules extension])
fi

AX_F90_LIBRARY_SETUP(ppm,[],ppm_module_finalize.$ax_cv_f90_modext,libppm*,[
    use ppm_module_finalize
    integer :: info
    call ppm_finalize(info)
   ])

#have_ppm=no
#AC_SEARCH_LIBS([ppm_module_finalize_mp_ppm_finalize],[ppm],[have_ppm=yes])
#if test "x${have_ppm}" = xno; then
#   echo "------------------------------------------"
#   echo " The PPM library is required to build numerics. "
#   echo " Either the library is missing or incorrectly compiled/linked."
#   echo "------------------------------------------"
#   (exit 1); exit 1;
#fi  

#AC_CHECK_LIB([metis],[METIS_PartGraphRecursive])
#have_metis=no
#AC_SEARCH_LIBS([METIS_PartGraphRecursive],[metis],[have_metis=yes])
#if test "x${have_metis}" = xno; then
#  echo "------------------------------------------"
#  echo " The METIS library is required to build PPM. "
#  echo " Either the library is missing or incorrectly compiled/linked."
#  echo " Stopping...   Please download at:"
#  echo " http://glaros.dtc.umn.edu/gkhome/metis/metis/download"
#  echo "------------------------------------------"
#  (exit 1); exit 1;
#else
#  echo "#define __METIS" >> ./src/ppm_define.h
#fi

use_mpi=no
AC_ARG_ENABLE([mpi],,[use_mpi=yes])
if test "x${use_mpi}" = xyes; then
  AC_CHECK_PROG([FC],[mpif90],[mpif90])
  AC_CHECK_LIB([mpi],[MPI_Init])
  have_mpi=no
  AC_SEARCH_LIBS([MPI_Init],[mpi],[have_mpi=yes])
  if test "x${have_mpi}" = xno; then
    echo "------------------------------------------"
    echo " The MPI library was not found."
    echo " Creating Makefile without parallel support."
    echo " Please download at:"
    echo " http://www.open-mpi.org/software/ "
    echo "------------------------------------------"
  else
    echo "#define __MPI" >> ./src/ppm_define.h
  fi
else
  echo "MPI was not enabled. Building non-parallel version of PPM." 
fi

# Check to see if fftw3 was enabled with -enable-fftw3
#use_fftw3=no
#AC_ARG_ENABLE([fftw3],,[use_fftw3=yes])

# If it was enabled than check that the library exists
#if test "x${use_fftw3}" = xyes; then
#  have_fftw=no
#  AC_CHECK_LIB([fftw3],[dfftw_execute])
#  AC_SEARCH_LIBS([dfftw_execute],[fftw3],[have_fftw=yes])
#  if test "x${have_fftw}" = xno; then
#    echo "------------------------------------------"
#    echo " You are trying to use FFTW to build PPM. "
#    echo " Either the library is missing or incorrectly compiled/linked."
#    echo " It can be downloaded at: http://www.fftw.org/download.html"
#    echo "------------------------------------------"
#    (exit 1); exit 1;
#  fi
#  echo "#define __FFTW" >> ./src/ppm_define.h
#fi

# Check to see if linux was enabled with -enable-linux
use_linux=no
AC_ARG_ENABLE([linux],,[use_linux=yes])
if test "x${use_linux}" = xyes; then
  echo "#define __Linux" >> ./src/ppm_define.h
fi

# Check to see if etime was enabled with -enable-etime
use_etime=no
AC_ARG_ENABLE([etime],,[use_etime=yes])
if test "x${use_etime}" = xyes; then
  echo "#define __ETIME" >> ./src/ppm_define.h
fi

# Check to see if vector was enabled with -enable-vector
use_vector=no
AC_ARG_ENABLE([vector],,[use_vector=yes])
if test "x${use_vector}" = xyes; then
  echo "#define __VECTOR" >> ./src/ppm_define.h
fi

use_mathkeisan=no
AC_ARG_ENABLE([mathkeisan],,[use_mathkeisan=yes])
if test "x${use_mathkeisan}" = xyes; then
  echo "#define __MATHKEISAN" >> ./src/ppm_define.h
fi

use_sxf90=no
AC_ARG_ENABLE([sxf90],,[use_sxf90=yes])
if test "x${use_sxf90}" = xyes; then
  echo "#define __SXF90" >> ./src/ppm_define.h
fi

use_crayfishpack=no
AC_ARG_ENABLE([crayfishpack],,[use_crayfishpack=yes])
if test "x${use_crayfishpack}" = xyes; then
  echo "#define __CRAYFISHPACK" >> ./src/ppm_define.h
fi

use_hypre=no
AC_ARG_ENABLE([hypre],,[use_hypre=yes])
if test "x${use_hypre}" = xyes; then
  echo "#define __HYPRE" >> ./src/ppm_define.h
fi

use_nomicroinstructions=no
AC_ARG_ENABLE([nomicroinstructions],,[use_nomicroinstructions=yes])
if test "x${use_nomicroinstructions}" = xyes; then
  echo "#define __NOMICROINSTRUCTIONS" >> ./src/ppm_define.h
fi

use_debug=no
AC_ARG_ENABLE([debug],,[use_debug=yes])
if test "x${use_debug}" = xyes; then
  AC_SUBST(DEBUG,-O0)
fi

# Switch language back to C for checking libraries
AC_LANG([C])

# Checks for header files
AC_HEADER_STDC
AC_CHECK_HEADERS([stdlib.h])

# Checks for typedefs and structures
AC_C_CONST
AC_C_INLINE

# Generate Makefiles
AC_CONFIG_FILES([Makefile])
AC_OUTPUT
