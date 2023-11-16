#!/usr/bin/bash

COMPILE_MODE     = release
USE_PROFILE      = 0
PAR_ND_LIBRARY   = ptscotch
SEQ_ND_LIBRARY   = scotch
USE_SYMPACK      = 1

 

# Different compiling and linking options.

# Edison is a CRAY XC30 machine

SUFFIX       = serverSYMPACK

 

# Compiler and tools

################################################################

CC           = /usr/bin/mpicc
CXX          = /usr/bin/mpicxx
FC           = /usr/bin/mpifort
LOADER       = /usr/bin/mpicxx

 

AR           = ar
ARFLAGS      = rvcu

# For System V based machine without ranlib, like Cray and SGI,

# use touch instead.

#RANLIB      = touch

RANLIB       = ranlib
 

CP           = cp
RM           = rm
RMFLAGS      = -f

################################################################

ENABLE_METIS=ON
ENABLE_PARMETIS=ON
METIS_LIB = ${METIS_DIR}/lib/libmetis.a
PARMETIS_LIB = ${METIS_DIR}/lib/libparmetis.a
METIS_INCLUDE = -I${METIS_DIR}/include
PARMETIS_INCLUDE = -I${METIS_DIR}/include
USE_METIS=1
USE_PARMETIS=1
USE_SCOTCH=1

# PEXSI directory
PEXSI_BUILD_DIR = ${PEXSI_DIR}/build


# Required libraries directories

# Includes
PEXSI_INCLUDE    = -I${PEXSI_DIR}/include
DSUPERLU_INCLUDE = -I${SUPERLU_DIR}/include
INCLUDES         = ${PEXSI_INCLUDE} ${DSUPERLU_INCLUDE}


# Libraries
DSUPERLU_LIB     = ${SUPERLU_DIR}/lib/libsuperlu_dist.a
PEXSI_LIB        = ${PEXSI_DIR}/src/libpexsi_${SUFFIX}.a
MKL_LIB          = -mkl


# Graph partitioning libraries

METIS_LIB        = -L${METIS_DIR} -lmetis
PARMETIS_LIB     = -L${PARMETIS_DIR} -lparmetis
SCOTCH_LIB       = -L${SCOTCH_DIR} -lscotchmetis -lscotch -lscotcherr
PTSCOTCH_LIB     = -L${PTSCOTCH_DIR} -lptscotchparmetis -lptscotch -lptscotcherr -lscotch


LIBS  = ${PEXSI_LIB} ${DSUPERLU_LIB} ${PAR_ND_LIB} ${SEQ_ND_LIB} ${MKL_LIB} ${PARMETIS_LIB} ${METIS_LIB} ${SCOTCH_LIB} ${PTSCOTCH_LIB}

COMPILE_DEF  += -DAdd_

#CPPFLAG = -std=c++11

 

#symPACK related definitions

include ${SYMPACK_DIR_2}/include/sympack.mak

CPPFLAG += ${SYMPACK_INCLUDE_2} ${SYMPACK_INCLUDE}
LIBS += ${METIS_LIB} ${PARMETIS_LIB}
CPPFLAG += ${METIS_INCLUDE} ${PARMETIS_INCLUDE}

LIBS+= ${SYMPACK_LIB_2} ${LAPACK_LIB} ${BLAS_LIB} ${GFORTRAN_LIB}

COMPILE_DEF  += -DWITH_SYMPACK

#CPPFLAG += -std=c++11

CFLAGS       = ${COMPILE_FLAG} ${PROFILE_FLAG} ${INCLUDES} #-std=c99
FFLAGS       = ${COMPILE_FLAG} ${PROFILE_FLAG} ${INCLUDES}
CXXFLAGS     = ${COMPILE_FLAG} ${CPPFLAG} ${PROFILE_FLAG} ${INCLUDES}
CCDEFS       = ${COMPILE_DEF}
CPPDEFS      = ${COMPILE_DEF}
LOADOPTS     = ${PROFILE_FLAG} ${LIBS} -Wl,--allow-multiple-definition
FLOADOPTS    = ${LIBS} -lstdc++ -Wl,--allow-multiple-definition

CCDEFS   += ${SYMPACK_DEFINITIONS}
CPPDEFS  += ${SYMPACK_DEFINITIONS}

