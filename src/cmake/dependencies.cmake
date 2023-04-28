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

# TIRPC
# ----
find_package(TIRPC REQUIRED)
include_directories(${TIRPC_INCLUDE_DIRS})

set(CMAKE_REQUIRED_INCLUDES ${TIRPC_INCLUDE_DIRS})
set(CMAKE_REQUIRED_LIBRARIES ${TIRPC_LIBRARY})
CHECK_CXX_SOURCE_COMPILES("
#include <stdint.h>
#include <rpc/types.h>
#include <rpc/xdr.h>
int main(int count, char** v){
	char buffer[15] = \"aaaaaaaaaaaaa\";
	XDR xdr;
	xdrmem_create(&xdr, buffer, 32, XDR_ENCODE);
	uint16_t a;
	uint32_t b;
	uint64_t c;
	xdr_uint16_t(&xdr, &a);
	xdr_uint32_t(&xdr, &b);
	xdr_uint64_t(&xdr, &c);
	return b;
}" HAVE_XDRUINTXX_T)

# without the standard names for the xdr functions, create aliases for the existing ones
if(NOT HAVE_XDRUINTXX_T)
	add_definitions(-Dxdr_uint16_t=xdr_u_int16_t)
	add_definitions(-Dxdr_uint32_t=xdr_u_int32_t)
	add_definitions(-Dxdr_uint64_t=xdr_u_int64_t)
endif()
