# check for header file
find_path(PARMETIS_INCLUDE_DIRS parmetis.h
	HINTS ${PARMETIS_DIR}/include $ENV{PARMETIS_DIR}/include
	PATH_SUFFIXES parmetis
	DOC "directory where the PARMETIS header is located")

# check for metis
find_library(METIS_LIBRARY
	NAMES metis
	HINTS ${METIS_DIR}/lib $ENV{METIS_DIR}/lib
	NO_DEFAULT_PATH
	DOC "the METIS library")
find_library(METIS_LIBRARY
	NAMES metis
	DOC "the METIS library")

# check for parmetis
find_library(PARMETIS_LIBRARY
	NAMES parmetis
	HINTS ${PARMETIS_DIR}/lib $ENV{PARMETIS_DIR}/lib
	NO_DEFAULT_PATH
	DOC "the PARMETIS library")
find_library(PARMETIS_LIBRARY
	NAMES parmetis
	DOC "the PARMETIS library")

set(PARMETIS_LIBRARIES ${PARMETIS_LIBRARY} ${METIS_LIBRARY})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PARMETIS
	PARMETIS_INCLUDE_DIR PARMETIS_LIBRARIES)
