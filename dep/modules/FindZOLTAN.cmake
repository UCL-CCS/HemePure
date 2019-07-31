# check for header file
find_path(ZOLTAN_INCLUDE_DIRS zoltan_cpp.h
	HINTS ${ZOLTAN_DIR}/include $ENV{ZOLTAN_DIR}/include
	PATH_SUFFIXES zoltan
	DOC "directory where the ZOLTAN header is located")

# check for zoltan
find_library(ZOLTAN_LIBRARY
	NAMES zoltan
	HINTS ${ZOLTAN_DIR}/lib $ENV{ZOLTAN_DIR}/lib
	NO_DEFAULT_PATH
	DOC "the ZOLTAN library")
find_library(ZOLTAN_LIBRARY
	NAMES zoltan
	DOC "the ZOLTAN library")

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(ZOLTAN
	ZOLTAN_INCLUDE_DIR ZOLTAN_LIBRARY)
