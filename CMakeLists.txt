#---------------------------------------------------
# OpenMM XTB Plugin
#----------------------------------------------------

CMAKE_MINIMUM_REQUIRED(VERSION 3.14)

# We need to know where OpenMM is installed so we can access the headers and libraries.
SET(OPENMM_DIR "${OPENMM_ROOT}" CACHE PATH "Where OpenMM is installed")
INCLUDE_DIRECTORIES("${OPENMM_DIR}/include")
message("OpenMM: ${OPENMM_DIR}")
LINK_DIRECTORIES("${OPENMM_DIR}/lib" "${OPENMM_DIR}/lib/plugins")

# We also need to know where XTB is installed.
FIND_FILE(XTB_HEADER_PATH xtb.h
        PATHS
        "${XTB_ROOT}/include"
        "${XTB_ROOT}/include/xtb"
        NO_DEFAULT_PATH
)
FIND_LIBRARY(XTB_LIBRARY_PATH xtb
        PATHS
        "${XTB_ROOT}/lib"
        "${XTB_ROOT}/lib/xtb"
        NO_DEFAULT_PATH
)
CMAKE_PATH(GET XTB_LIBRARY_PATH PARENT_PATH XTB_LIBRARY_DIR)
CMAKE_PATH(GET XTB_HEADER_PATH PARENT_PATH XTB_INCLUDE_DIR)
INCLUDE_DIRECTORIES("${XTB_INCLUDE_DIR}")
LINK_DIRECTORIES("${XTB_LIBRARY_DIR}")

# Specify the C++ version we are building for.
SET (CMAKE_CXX_STANDARD 20)

# Set flags for linking on mac
IF(APPLE)
    SET(CMAKE_INSTALL_NAME_DIR "@rpath")
    SET(EXTRA_COMPILE_FLAGS "-stdlib=libc++")
ENDIF(APPLE)

# Select where to install
IF(${CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT})
    IF(WIN32)
        SET(CMAKE_INSTALL_PREFIX "$ENV{ProgramFiles}/OpenMM" CACHE PATH "Where to install the plugin" FORCE)
    ELSE(WIN32)
        SET(CMAKE_INSTALL_PREFIX "${INSTALL_ROOT}" CACHE PATH "Where to install the plugin" FORCE)
    ENDIF(WIN32)
ENDIF(${CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT})


# The source is organized into subdirectories, but we handle them all from
# this CMakeLists file rather than letting CMake visit them as SUBDIRS.
SET(XTB_PLUGIN_SOURCE_SUBDIRS openmmapi serialization)

# Set the library name
SET(XTB_LIBRARY_NAME OpenMMXTB)
SET(SHARED_XTB_TARGET ${XTB_LIBRARY_NAME})

# These are all the places to search for header files which are to be part of the API.
SET(API_INCLUDE_DIRS "openmmapi/include" "openmmapi/include/internal")

# Locate header files.
SET(API_INCLUDE_FILES)
FOREACH(dir ${API_INCLUDE_DIRS})
    FILE(GLOB fullpaths ${dir}/*.h)
    SET(API_INCLUDE_FILES ${API_INCLUDE_FILES} ${fullpaths})
ENDFOREACH(dir)

# Collect up source files
SET(SOURCE_FILES) # empty
SET(SOURCE_INCLUDE_FILES)
FOREACH(subdir ${XTB_PLUGIN_SOURCE_SUBDIRS})
    FILE(GLOB src_files  ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/src/*.cpp)
    FILE(GLOB incl_files ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/src/*.h)
    SET(SOURCE_FILES         ${SOURCE_FILES}         ${src_files})   #append
    SET(SOURCE_INCLUDE_FILES ${SOURCE_INCLUDE_FILES} ${incl_files})

    ## Make sure we find these locally before looking in OpenMM/include if
    ## OpenMM was previously installed there.
    INCLUDE_DIRECTORIES(BEFORE ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/include)
ENDFOREACH(subdir)

# Create the library.

ADD_LIBRARY(${SHARED_XTB_TARGET} SHARED ${SOURCE_FILES} ${SOURCE_INCLUDE_FILES} ${API_INCLUDE_FILES})
SET_TARGET_PROPERTIES(${SHARED_XTB_TARGET}
    PROPERTIES COMPILE_FLAGS "-DXTB_BUILDING_SHARED_LIBRARY ${EXTRA_COMPILE_FLAGS}"
    LINK_FLAGS "${EXTRA_COMPILE_FLAGS}")
TARGET_LINK_LIBRARIES(${SHARED_XTB_TARGET} OpenMM xtb)
INSTALL_TARGETS(/lib RUNTIME_DIRECTORY /lib ${SHARED_XTB_TARGET})
TARGET_COMPILE_OPTIONS(${SHARED_XTB_TARGET} PUBLIC -pg)

# install headers
FILE(GLOB API_ONLY_INCLUDE_FILES "openmmapi/include/*.h")
INSTALL (FILES ${API_ONLY_INCLUDE_FILES} DESTINATION include)
FILE(GLOB API_ONLY_INCLUDE_FILES_INTERNAL "openmmapi/include/internal/*.h")
INSTALL (FILES ${API_ONLY_INCLUDE_FILES_INTERNAL} DESTINATION include/internal)

# Enable testing

ENABLE_TESTING()
ADD_SUBDIRECTORY(tests)
ADD_SUBDIRECTORY(serialization/tests)

# Build the Python API

#FIND_PROGRAM(PYTHON_EXECUTABLE python)
#FIND_PROGRAM(PIP_EXECUTABLE pip)
FIND_PROGRAM(PYTHON_EXECUTABLE python
        PATHS
        "${PYTHON_ROOT}/bin"
        "${PYTHON_ROOT}/envs/*/bin"
        NO_DEFAULT_PATH)

FIND_PROGRAM(PIP_EXECUTABLE pip
        PATHS
        "${PYTHON_ROOT}/bin"
        "${PYTHON_ROOT}/envs/*/bin"
        NO_DEFAULT_PATH)

FIND_PROGRAM(SWIG_EXECUTABLE swig)
IF(PYTHON_EXECUTABLE AND SWIG_EXECUTABLE)
    SET(XTB_BUILD_PYTHON_WRAPPERS ON CACHE BOOL "Build wrappers for Python")
ELSE(PYTHON_EXECUTABLE AND SWIG_EXECUTABLE)
    SET(XTB_BUILD_PYTHON_WRAPPERS OFF CACHE BOOL "Build wrappers for Python")
ENDIF(PYTHON_EXECUTABLE AND SWIG_EXECUTABLE)
IF(XTB_BUILD_PYTHON_WRAPPERS)
    ADD_SUBDIRECTORY(python)
ENDIF(XTB_BUILD_PYTHON_WRAPPERS)
