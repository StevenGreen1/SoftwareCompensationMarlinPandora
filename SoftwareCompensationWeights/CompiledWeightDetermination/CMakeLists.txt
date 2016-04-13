########################################################
# cmake file for building CalorimeterAnalysis
# @author Steve Green
CMAKE_MINIMUM_REQUIRED(VERSION 2.6 FATAL_ERROR)
########################################################

# project name
PROJECT( SoftCompWeightDetermination )

# project version
SET( ${PROJECT_NAME}_VERSION_MAJOR 0 )
SET( ${PROJECT_NAME}_VERSION_MINOR 0 )
SET( ${PROJECT_NAME}_VERSION_PATCH 0 )

### DEPENDENCIES ############################################################

FIND_PACKAGE( ILCUTIL COMPONENTS ILCSOFT_CMAKE_MODULES REQUIRED )

# load default settings from ILCSOFT_CMAKE_MODULES
INCLUDE( ilcsoft_default_settings )

FIND_PACKAGE( ROOT REQUIRED )

FOREACH( pkg ROOT )
    IF( ${pkg}_FOUND )
        INCLUDE_DIRECTORIES( ${${pkg}_INCLUDE_DIRS} )
        LINK_LIBRARIES( ${${pkg}_LIBRARIES} /afs/desy.de/project/ilcsoft/sw/x86_64_gcc44_sl6/v01-17-05/root/5.34.10/lib/libMinuit.so )
        ADD_DEFINITIONS ( ${${pkg}_DEFINITIONS} )
    ENDIF()
ENDFOREACH()

INCLUDE_DIRECTORIES(BEFORE ./include ) # use BEFORE to fix problem on mac osx (case insensitive filesystem)

### LIBRARY AND TOOLS #######################################################

ADD_EXECUTABLE(SoftCompWeightDetermination ./src/SoftCompWeightDetermination.cxx)

INSTALL( TARGETS CalorimeterAnalysis DESTINATION bin PERMISSIONS
    OWNER_READ OWNER_WRITE OWNER_EXECUTE
    GROUP_READ GROUP_EXECUTE
    WORLD_READ WORLD_EXECUTE )
