# check for header file
find_path(CTEMPLATE_INCLUDE_DIRS ctemplate/template.h
	HINTS ${CTEMPLATE_DIR}/include $ENV{CTEMPLATE_DIR}/include
	PATH_SUFFIXES ctemplate
	DOC "directory where the CTEMPLATE header is located")

# check for ctemplate
find_library(CTEMPLATE_LIBRARY
	NAMES ctemplate
	HINTS ${CTEMPLATE_DIR}/lib $ENV{CTEMPLATE_DIR}/lib
	NO_DEFAULT_PATH
	DOC "the CTEMPLATE library")
find_library(CTEMPLATE_LIBRARY
	NAMES ctemplate
	DOC "the CTEMPLATE library")

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(CTEMPLATE
	CTEMPLATE_INCLUDE_DIR CTEMPLATE_LIBRARY)
