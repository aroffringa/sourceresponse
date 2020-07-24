cmake_minimum_required(VERSION 2.8.8)

project(sourceresponse)

add_definitions(-DAOPROJECT)

if(PORTABLE)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -ggdb -Wvla -Wall -DNDEBUG -std=c++11")
else()
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -ggdb -Wvla -Wall -Wextra -DNDEBUG -march=native -std=c++11")
endif(PORTABLE)

# Casacore has a separate CMake file in this directory
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMake)

find_package(HDF5 COMPONENTS C CXX REQUIRED)

SET(CASACORE_MAKE_REQUIRED_EXTERNALS_OPTIONAL TRUE)
find_package(Casacore REQUIRED COMPONENTS casa ms tables measures fits)

find_package(CFITSIO REQUIRED)

#Prevent accidentally finding old BoostConfig.cmake file from casapy
set(Boost_NO_BOOST_CMAKE ON)
find_package(Boost COMPONENTS date_time system filesystem program_options REQUIRED)
find_library(PTHREAD_LIB pthread REQUIRED)
find_library(DL_LIB dl REQUIRED)
find_library(GSL_LIB NAMES gsl)
find_path(GSL_INCLUDE_DIR NAMES gsl/gsl_version.h)
find_library(GSL_CBLAS_LIB NAMES gslcblas)

include_directories(${CMAKE_SOURCE_DIR}/external)

find_path(LOFAR_STATION_RESPONSE_DIR NAMES StationResponse/ITRFConverter.h)
find_library(LOFAR_STATION_RESPONSE_LIB stationresponse)
if(LOFAR_STATION_RESPONSE_DIR AND LOFAR_STATION_RESPONSE_LIB)
	include_directories(${LOFAR_STATION_RESPONSE_DIR})
  set(LBEAM_LIBS ${LOFAR_STATION_RESPONSE_LIB})
  set(LBEAM_FILES lofar/lbeamevaluator.cpp lofar/lmspredicter.cpp aterms/lofarbeamterm.cpp)
  set(HAVE_LOFAR_BEAM True)
  add_definitions(-DHAVE_LOFAR_BEAM)
  message(STATUS "LOFAR beam library found.")
else()
  message(FATAL_ERROR "LOFAR beam library not found.")
endif(LOFAR_STATION_RESPONSE_DIR AND LOFAR_STATION_RESPONSE_LIB)

include_directories(${CASACORE_INCLUDE_DIRS})
include_directories(${Boost_INCLUDE_DIR})
include_directories(${CFITSIO_INCLUDE_DIR})
include_directories(${FFTW3_INCLUDE_DIR})
include_directories(${HDF5_INCLUDE_DIRS})

if(GSL_LIB AND GSL_CBLAS_LIB AND GSL_INCLUDE_DIR)
  include_directories(${GSL_INCLUDE_DIR})
endif(GSL_LIB AND GSL_CBLAS_LIB AND GSL_INCLUDE_DIR)

# use, i.e. don't skip the full RPATH for the build tree
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)
# when building, don't use the install RPATH already
# (but later on when installing)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE) 
SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
# the RPATH to be used when installing, but only if it's not a system directory
LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
IF("${isSystemDir}" STREQUAL "-1")
   SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
ENDIF("${isSystemDir}" STREQUAL "-1")

add_executable(sourceresponse sourceresponse.cpp model/model.cpp nlplfitter.cpp polynomialfitter.cpp)
target_link_libraries(sourceresponse ${CFITSIO_LIBRARY} ${CASACORE_LIBRARIES} ${GSL_LIB} ${GSL_CBLAS_LIB} ${Boost_SYSTEM_LIBRARY} ${Boost_DATE_TIME_LIBRARY} ${LBEAM_LIBS})

if(NOT GSL_CBLAS_LIB)
  message(WARNING "GSL CBLAS lib was not found. GSL needs CBLAS: disabling GSL.")
endif(NOT GSL_CBLAS_LIB)

if(GSL_LIB AND GSL_CBLAS_LIB AND GSL_INCLUDE_DIR)
  add_definitions(-DHAVE_GSL)
#  target_link_libraries(.. ${GSL_LIB} ${GSL_CBLAS_LIB})
else(GSL_LIB AND GSL_CBLAS_LIB AND GSL_INCLUDE_DIR)
  message(WARNING "GNU scientific library (gsl) not found: some features will be disabled (this is not an error).")
endif(GSL_LIB AND GSL_CBLAS_LIB AND GSL_INCLUDE_DIR)

message(STATUS "Flags passed to C++ compiler: " ${CMAKE_CXX_FLAGS})
