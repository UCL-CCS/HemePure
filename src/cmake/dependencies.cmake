## PTSCOTCH
## --------
#find_package(PTSCOTCH REQUIRED)
#include_directories(${PTSCOTCH_INCLUDE_DIRS})

# PARMETIS
# --------
find_package(PARMETIS REQUIRED)
include_directories(${PARMETIS_INCLUDE_DIRS})

## ZOLTAN
## ------
#find_package(ZOLTAN REQUIRED)
#include_directories(${ZOLTAN_INCLUDE_DIRS})

# TINYXML
# -------
find_package(TINYXML REQUIRED)
option(TIXML_USE_STL "Use STL with TIXML" ON)
if(TIXML_USE_STL)
	add_definitions(-DTIXML_USE_STL)
endif()
include_directories(${TINYXML_INCLUDE_DIR})

# CTEMPLATE
# ---------
find_package(CTEMPLATE REQUIRED)
include_directories(${CTEMPLATE_INCLUDE_DIR})

# BOOST
# -----
find_package(Boost 1.54 REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

# ZLIB
# ----
find_package(ZLIB REQUIRED)
include_directories(${ZLIB_INCLUDE_DIR})
