##################################################
# LATTICE DMRG PROJECT
##################################################
CMAKE_MINIMUM_REQUIRED(VERSION 3.10.0)
PROJECT(dmrg Fortran)
SET(VERSION 4.1.0)


MESSAGE(STATUS "OS: ${CMAKE_SYSTEM_NAME} ${CMAKE_SYSTEM_VERSION}")
SET(BUILD_TYPE "RELEASE" CACHE STRING "Build type")
STRING(TOUPPER "${BUILD_TYPE}" BUILD_TYPE)
MESSAGE(STATUS "BUILD_TYPE=${BUILD_TYPE}")
SET(CMAKE_BUILD_TYPE ${BUILD_TYPE} CACHE INTERNAL "Build type" FORCE)
#
SET(PRECISION "DBLE" CACHE STRING "Variable Precision type")
STRING(TOUPPER "${PRECISION}" PRECISION)
MESSAGE(STATUS "PRECISION=${PRECISION}")
SET(CMAKE_PRECISION ${PRECISION} CACHE INTERNAL "Variable Precision type" FORCE)


# Add our local CMAKE modules to the module path
SET(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
INCLUDE(${CMAKE_MODULE_PATH}/ColorsMsg.cmake)

# Set some options the user may choose
OPTION(USE_MPI "Set MPI pre-compilation flag" ON)
OPTION(PROFILE "Set Code Profiling" ON)
OPTION(API_EXPORT "Export full interface" OFF)
OPTION(VERBOSE "Build a verbose Makefile" OFF)
OPTION(WITH_BLAS_LAPACK "Use internal Blas/Lapack" OFF)
IF(VERBOSE)
  SET(CMAKE_VERBOSE_MAKEFILE ON)
ELSE(VERBOSE)
  SET(CMAKE_VERBOSE_MAKEFILE OFF)
ENDIF(VERBOSE)

MESSAGE(STATUS "USE_MPI=${USE_MPI}")
MESSAGE(STATUS "PROFILE=${PROFILE}")
MESSAGE(STATUS "API_EXPORT=${API_EXPORT}")

############################################################
# Look for MPI if needed
############################################################
INCLUDE(${CMAKE_MODULE_PATH}/MpiConfig.cmake)


############################################################
# Look for required libraries: scifor
# (includes Blas/Lapack/Scalapack/Quadpack/Minpack/FFTpack/Arpack/...) 
############################################################
# Find Scifor
FIND_PACKAGE(PkgConfig)
PKG_SEARCH_MODULE(SCIFOR REQUIRED scifor)
IF(${SCIFOR_FOUND})
  MESSAGE(STATUS "SciFortran found ${SCIFOR_INCLUDE_DIRS}")
ELSE()
  MESSAGE(FATAL_ERROR "
Can not find SciFortran in this system. Check SciFor installation. 
Try adding the provided PKG-CONFIG file to the environment:
$ export PKG_CONFIG_PATH=<path_to_scifor_dir>/etc:$PKG_CONFIG_PATH")
ENDIF()
INCLUDE_DIRECTORIES(BEFORE ${SCIFOR_INCLUDE_DIRS})

# FIND_PACKAGE(Gperftools OPTIONAL)
# if(${Gperftools_FOUND})
#   MESSAGE(STATUS "GperfTools found ${Gperftools_INCLUDE_DIRS}")
#   INCLUDE_DIRECTORIES(BEFORE ${Gperftools_INCLUDE_DIRS})
# endif()


############################################################
# Configure Fortran Compiler Main options
############################################################
INCLUDE(${CMAKE_MODULE_PATH}/FortranCompilerConfig.cmake)




############################################################
# Check for .git and extract the current SHA1 revision.
############################################################
INCLUDE(${CMAKE_MODULE_PATH}/GitConfig.cmake)




############################################################
# Define the actual files and folders that make up the build
############################################################
INCLUDE(${CMAKE_MODULE_PATH}/MainConfig.cmake)



############################################################
############################################################
############################################################
# DEFINE HERE THE DETAILS OF THE LIBRARY CONFIGURATION 

# Define some directories
SET(LIB_SRC ${CMAKE_SOURCE_DIR}/src)
SET(LIB_ETC ${CMAKE_SOURCE_DIR}/etc)
SET(LIB_TEST ${CMAKE_SOURCE_DIR}/test)
SET(LIB_ENV ${CMAKE_SOURCE_DIR}/etc/environment_modules)
SET(LIB_TMP_INC ${CMAKE_BINARY_DIR}/include)
SET(LIB_TMP_ETC ${CMAKE_BINARY_DIR}/etc)
SET(LIB_TMP_VER ${CMAKE_BINARY_DIR}/version)
SET(LIB_TARGET_LIB ${CMAKE_INSTALL_PREFIX}/lib)
SET(LIB_TARGET_INC ${CMAKE_INSTALL_PREFIX}/include)
SET(LIB_TARGET_ETC ${CMAKE_INSTALL_PREFIX}/etc)
SET(LIB_TARGET_BIN ${CMAKE_INSTALL_PREFIX}/bin)
SET(LIB_TARGET_DIR ${CMAKE_INSTALL_PREFIX})

# Have the .mod files placed in the INC folder
SET(CMAKE_Fortran_MODULE_DIRECTORY ${LIB_TMP_INC})
INCLUDE_DIRECTORIES(BEFORE ${CMAKE_Fortran_MODULE_DIRECTORY})


ADD_SUBDIRECTORY(${LIB_SRC}/MATRIX)

ADD_SUBDIRECTORY(${LIB_SRC}/LIST)

ADD_SUBDIRECTORY(${LIB_SRC}/FOCK)

ADD_SUBDIRECTORY(${LIB_SRC}/BLOCKS)

ADD_SUBDIRECTORY(${LIB_SRC})


ADD_DEPENDENCIES(MATRIX_LIB
  AUX_LIB
)

ADD_DEPENDENCIES(FOCK_LIB
  AUX_LIB
)


ADD_DEPENDENCIES(LIST_LIB
  AUX_LIB
  MATRIX_LIB
)

ADD_DEPENDENCIES(BLOCKS_LIB
  AUX_LIB
  MATRIX_LIB
  FOCK_LIB
  LIST_LIB
)


ADD_DEPENDENCIES(DMRG_LIB
  AUX_LIB
  FOCK_LIB
  MATRIX_LIB
  LIST_LIB
  BLOCKS_LIB
  )

ADD_LIBRARY(dmrg STATIC
  $<TARGET_OBJECTS:AUX_LIB>
  $<TARGET_OBJECTS:FOCK_LIB>
  $<TARGET_OBJECTS:MATRIX_LIB>
  $<TARGET_OBJECTS:LIST_LIB>
  $<TARGET_OBJECTS:BLOCKS_LIB>
  $<TARGET_OBJECTS:DMRG_LIB>
  )




SET(PERMISSION_777
  OWNER_WRITE OWNER_READ OWNER_EXECUTE
  GROUP_WRITE GROUP_READ GROUP_EXECUTE
  WORLD_WRITE WORLD_READ WORLD_EXECUTE)


#Build the PKG-CONFIG file
SET(TMP_PKCONFIG_FILE ${LIB_TMP_ETC}/${PROJECT_NAME}.pc)
CONFIGURE_FILE( ${LIB_ETC}/${PROJECT_NAME}.pc.in ${TMP_PKCONFIG_FILE} @ONLY)

############################################################
############################################################
############################################################


INCLUDE(${CMAKE_MODULE_PATH}/PostConfig.cmake)
