# check for header file
find_path(TINYXML_INCLUDE_DIRS tinyxml.h
	HINTS ${TINYXML_DIR}/include $ENV{TINYXML_DIR}/include
	PATH_SUFFIXES tinyxml
	DOC "directory where the TINYXML header is located")

# check for tinyxml
find_library(TINYXML_LIBRARY
	NAMES tinyxml
	HINTS ${TINYXML_DIR}/lib $ENV{TINYXML_DIR}/lib
	NO_DEFAULT_PATH
	DOC "the TINYXML library")
find_library(TINYXML_LIBRARY
	NAMES tinyxml
	DOC "the TINYXML library")

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(TINYXML
	TINYXML_INCLUDE_DIR TINYXML_LIBRARY)
