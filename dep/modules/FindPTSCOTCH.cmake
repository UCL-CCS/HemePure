# check for header file
find_path(PTSCOTCH_INCLUDE_DIRS ptscotch.h
	HINTS ${PTSCOTCH_DIR}/include $ENV{PTSCOTCH_DIR}/include
	PATH_SUFFIXES ptscotch
	DOC "directory where the PTSCOTCH header is located")

# check for scotch
find_library(SCOTCH_LIBRARY
	NAMES scotch
	HINTS ${SCOTCH_DIR}/lib $ENV{SCOTCH_DIR}/lib
	NO_DEFAULT_PATH
	DOC "the SCOTCH library")
find_library(SCOTCH_LIBRARY
	NAMES scotch
	DOC "the SCOTCH library")

# check for scotcherr
find_library(SCOTCHERR_LIBRARY
	NAMES scotcherr
	HINTS ${SCOTCH_DIR}/lib $ENV{SCOTCH_DIR}/lib
	NO_DEFAULT_PATH
	DOC "the SCOTCH-ERROR library")
find_library(SCOTCHERR_LIBRARY
	NAMES scotcherr
	DOC "the SCOTCH-ERROR library")

# check for ptscotch
find_library(PTSCOTCH_LIBRARY
	NAMES ptscotch
	HINTS ${SCOTCH_DIR}/lib $ENV{SCOTCH_DIR}/lib
	NO_DEFAULT_PATH
	DOC "the PTSCOTCH library")
find_library(PTSCOTCH_LIBRARY
	NAMES ptscotch
	DOC "the PTSCOTCH library")

# check for ptscotcherr
find_library(PTSCOTCHERR_LIBRARY
	NAMES ptscotcherr
	HINTS ${SCOTCH_DIR}/lib $ENV{SCOTCH_DIR}/lib
	NO_DEFAULT_PATH
	DOC "the PTSCOTCH-ERROR library")
find_library(PTSCOTCHERR_LIBRARY
	NAMES ptscotcherr
	DOC "the PTSCOTCH-ERROR library")

# check for ptscotchmetis
find_library(PTSCOTCHPARMETIS_LIBRARY
	NAMES ptscotchparmetis
	HINTS ${SCOTCH_DIR}/lib $ENV{SCOTCH_DIR}/lib
	NO_DEFAULT_PATH
	DOC "the PTSCOTCHPARMETIS library")
find_library(PTSCOTCHPARMETIS_LIBRARY
	NAMES ptscotchparmetis
	DOC "the PTSCOTCHPARMETIS library")

#set(PTSCOTCH_LIBRARIES ${PTSCOTCHPARMETIS_LIBRARY} ${PTSCOTCH_LIBRARY} ${PTSCOTCHERR_LIBRARY} ${SCOTCH_LIBRARY})
 set(PTSCOTCH_LIBRARIES ${PTSCOTCH_LIBRARY} ${PTSCOTCHERR_LIBRARY} ${SCOTCH_LIBRARY})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PTSCOTCH
	PTSCOTCH_INCLUDE_DIR PTSCOTCH_LIBRARIES)
