find_package(PkgConfig QUIET)
pkg_check_modules(PC_TIRPC libtirpc)

# check for header file
find_path(TIRPC_INCLUDE_DIRS
	NAMES netconfig.h
	HINTS ${TIRPC_DIR}/include $ENV{TIRPC_DIR}/include ${PC_TIRPC_INCLUDE_DIRS}
	PATH_SUFFIXES tirpc
	DOC "directory where the TIRPC header is located")

# check for tirpc
find_library(TIRPC_LIBRARY
	NAMES tirpc
	HINTS ${TIRPC_DIR}/lib $ENV{TIRPC_DIR}/lib ${PC_TIRPC_LIBRARY_DIRS}
	NO_DEFAULT_PATH
	DOC "the TIRPC library")
find_library(TIRPC_LIBRARY
	NAMES tirpc
	DOC "the TIRPC library")

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(TIRPC
	REQUIRED_VARS TIRPC_INCLUDE_DIRS TIRPC_LIBRARY
	VERSION_VAR TIRPC_VERSION)
