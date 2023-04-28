cmake_minimum_required(VERSION 3.1)

project(HemeLBDependencies)
include(ExternalProject)

set(HEMELB_DEPENDENCIES_PATH ${CMAKE_CURRENT_SOURCE_DIR})
set(HEMELB_DEPENDENCIES_INSTALL_PATH "${HEMELB_DEPENDENCIES_PATH}/install"
	CACHE STRING "path to install dependencies")
set(HEMELB_SUBPROJECT_MAKE_JOBS 1
	CACHE INTEGER "number of jobs to use for subproject build steps")
set(CMAKE_MODULE_PATH "${HEMELB_DEPENDENCIES_PATH}/modules/" ${CMAKE_MODULE_PATH})
set(CPPUNIT_CONFIGURE_OPTIONS ""
	CACHE STRING "extra configuration options for CPPUNIT")
option(CPPUNIT_PATCH_LDL "add -ldl option to cppunit" OFF)
option(CPPUNIT_PATCH_DYNAMIC "add -dynamic option to cppunit; only works if CPPUNIT_PATCH_LDL is OFF" OFF)
set(CTEMPLATE_CONFIGURE_OPTIONS ""
	CACHE STRING "extra configurations options for CTEMPLATE")

find_package(MPI REQUIRED)
set(CMAKE_CXX_COMPILE_FLAGS "${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS}")
set(CMAKE_CXX_LINK_FLAGS "${CMAKE_CXX_LINK_FLAGS} ${CMAKE_CXX_LINK_FLAGS}")
include_directories(${MPI_INCLUDE_PATH})

# TINYXML
# -------
message("Building TINYXML")

option(TIXML_USE_STL "Use STL with TIXML" ON)

if(TIXML_USE_STL)
	add_definitions(-DTIXML_USE_STL)
endif()

find_file(TINYXML_TARBALL tinyxml_2_6_2.tar.gz
	DOC "Path to download TINYXML (can be url http://)"
	PATHS ${HEMELB_DEPENDENCIES_PATH}/src)

if(NOT TINYXML_TARBALL)
	message("No TINYXML source found, will download")
	set(TINYXML_TARBALL http://sourceforge.net/projects/tinyxml/files/tinyxml/2.6.2/tinyxml_2_6_2.tar.gz
		CACHE STRING "Path to download TinyXML (can be local file://)" FORCE)
endif()

ExternalProject_Add(TinyXML
	INSTALL_DIR			${HEMELB_DEPENDENCIES_INSTALL_PATH}
	URL					${TINYXML_TARBALL}
	CONFIGURE_COMMAND	cmake <SOURCE_DIR>
						-DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
						-DTIXML_USE_STL=${TIXML_USE_STL}
						-DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
						-DCMAKE_CXX_FLAGS_RELEASE=${CMAKE_CXX_FLAGS_RELEASE}
	PATCH_COMMAND		cp ${HEMELB_DEPENDENCIES_PATH}/patches/tinyxml.cmake CMakeLists.txt
	BUILD_COMMAND		make -j${HEMELB_SUBPROJECT_MAKE_JOBS})

# ALL
# ---
message("Building ALL")

find_file(ALL_TARBALL ALL_current.tar.gz
	PATHS ${HEMELB_DEPENDENCIES_PATH}/src)

ExternalProject_Add(ALL
	INSTALL_DIR			${HEMELB_DEPENDENCIES_INSTALL_PATH}
	URL					${ALL_TARBALL}
	CONFIGURE_COMMAND	cmake <SOURCE_DIR>
						-DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
						-DCMAKE_C_COMPILER=${MPI_C_COMPILER}
						-DCMAKE_CXX_COMPILER=${MPI_CXX_COMPILER}
	BUILD_COMMAND		make -j${HEMELB_SUBPROJECT_MAKE_JOBS})

# PARMETIS
# --------
message("Building PARMETIS")

find_file(PARMETIS_TARBALL parmetis-4.0.2.tar.gz
	DOC "Path to download PARMETIS (can be url http://)"
	PATHS ${HEMELB_DEPENDENCIES_PATH}/src)

if(NOT PARMETIS_TARBALL)
	message("No PARMETIS source found, will download")
	set(PARMETIS_TARBALL http://glaros.dtc.umn.edu/gkhome/fetch/sw/parmetis/parmetis-4.0.2.tar.gz
		CACHE STRING "Path to download PARMETIS (can be local file://)" FORCE)
endif()

set(PARMETIS_CC  ${MPI_C_COMPILER}   CACHE STRING   "MPI_C compiler to use for PARMETIS, leave blank to let PARMETIS guess")
set(PARMETIS_CXX ${MPI_CXX_COMPILER} CACHE STRING "MPI_CXX compiler to use for PARMETIS, leave blank to let PARMETIS guess")

if(PARMETIS_CC)
	set(PARMETIS_CC_OPTION	cc=${PARMETIS_CC})
endif()

if(PARMETIS_CXX)
	set(PARMETIS_CXX_OPTION	cxx=${PARMETIS_CXX})
endif()

ExternalProject_Add(ParMETIS
	INSTALL_DIR			${HEMELB_DEPENDENCIES_INSTALL_PATH}
	URL					${PARMETIS_TARBALL}
	CONFIGURE_COMMAND	make config prefix=<INSTALL_DIR> ${PARMETIS_CC_OPTION} ${PARMETIS_CXX_OPTION} && cd metis &&
						make config prefix=<INSTALL_DIR> ${PARMETIS_CC_OPTION} ${PARMETIS_CXX_OPTION}
	BUILD_COMMAND		make -j${HEMELB_SUBPROJECT_MAKE_JOBS} MAKEFLAGS= # workaround for bug in ParMETIS makefile treating switches as targets
	INSTALL_COMMAND		make install MAKEFLAGS= && cd metis && make install MAKEFLAGS=
	BUILD_IN_SOURCE		1)

## SCOTCH
## ------
#message("Building SCOTCH")
#
#find_file(SCOTCH_TARBALL scotch_6.0.4.tar.gz
#	DOC "Path to download SCOTCH (can be url http://)"
#	PATHS ${HEMELB_DEPENDENCIES_PATH}/src)
#
#if(NOT SCOTCH_TARBALL)
#	message("No SCOTCH source found, will download")
#	set(SCOTCH_TARBALL http://gforge.inria.fr/frs/download.php/file/34618/scotch_6.0.4.tar.gz
#		CACHE STRING "Path to download SCOTCH (can be local file://)" FORCE)
#endif()
#
#ExternalProject_Add(scotch
#	INSTALL_DIR			${HEMELB_DEPENDENCIES_INSTALL_PATH}
#	URL					${SCOTCH_TARBALL}
#	CONFIGURE_COMMAND	cp ${HEMELB_DEPENDENCIES_PATH}/stuff/Makefile.inc <SOURCE_DIR>/src/Makefile.inc &&
#						cp ${HEMELB_DEPENDENCIES_PATH}/stuff/parmetis_dgraph_part.c <SOURCE_DIR>/src/libscotchmetis/
#	BUILD_COMMAND		cd src && make ptscotch -j${HEMELB_SUBPROJECT_MAKE_JOBS}
#	INSTALL_COMMAND		cd src && make prefix=<INSTALL_DIR> install
#	BUILD_IN_SOURCE		1)

## ZOLTAN
## --------
#message("Building ZOLTAN")
#
#find_file(ZOLTAN_TARBALL zoltan_distrib_v3.83.tar.gz
#	DOC "Path to download ZOLTAN (can be url http://)"
#	PATHS ${HEMELB_DEPENDENCIES_PATH}/src)
#
#if(NOT ZOLTAN_TARBALL)
#	message("No ZOLTAN source found, will download")
#	set(PARMETIS_TARBALL http://www.cs.sandia.gov/~kddevin/Zoltan_Distributions/zoltan_distrib_v3.83.tar.gz
#		CACHE STRING "Path to download ZOLTAN (can be local file://)" FORCE)
#endif()
#
#ExternalProject_Add(zoltan
#	INSTALL_DIR			${HEMELB_DEPENDENCIES_INSTALL_PATH}
#	URL					${ZOLTAN_TARBALL}
#	CONFIGURE_COMMAND	<SOURCE_DIR>/configure --prefix=<INSTALL_DIR>
#						--enable-mpi --with-mpi-compilers --with-gnumake
##						--with-scotch
##						--with-scotch-incdir="${HEMELB_DEPENDENCIES_INSTALL_PATH}/include"
##						--with-scotch-libdir="${HEMELB_DEPENDENCIES_INSTALL_PATH}/lib"
#						--with-parmetis
#						--with-parmetis-incdir="${HEMELB_DEPENDENCIES_INSTALL_PATH}/include"
#						--with-parmetis-libdir="${HEMELB_DEPENDENCIES_INSTALL_PATH}/lib"
#						--with-id-type=uint
#	BUILD_COMMAND		make everything -j${HEMELB_SUBPROJECT_MAKE_JOBS}
#	INSTALL_COMMAND		make install
#	BUILD_IN_SOURCE 0)

## CPPUNIT
## -------
#message("Building CPPUNIT")
#
#find_file(CPPUNIT_TARBALL cppunit-1.12.1.tar.gz
#	DOC "Path to download CPPUNIT (can be url http://)"
#	PATHS ${HEMELB_DEPENDENCIES_PATH}/src)
#
#if(NOT CPPUNIT_TARBALL)
#	message("No CPPUNIT source found, will download")
#	set(CPPUNIT_TARBALL http://downloads.sourceforge.net/cppunit/cppunit-1.12.1.tar.gz
#		CACHE STRING "Path to download CPPUNIT (can be local file://)" FORCE)
#endif()
#
#set(RECONFIGURE_CPPUNIT_DEFAULT OFF)
#option(RECONFIGURE_CPPUNIT "Include Reconfigure step for CPPUNIT" ${RECONFIGURE_CPPUNIT_DEFAULT})
#
#if (CPPUNIT_PATCH_LDL)
#	set(PATCH_COMMAND_LDL patch -p1 < ${HEMELB_DEPENDENCIES_PATH}/patches/cppunit_ldl.diff)
#elseif (CPPUNIT_PATCH_DYNAMIC)
#	set(PATCH_COMMAND_LDL patch -p1 < ${HEMELB_DEPENDENCIES_PATH}/patches/cppunit_dynamic.diff)
#else()
#	set(PATCH_COMMAND_LDL echo noldl)
#endif()
#
#if (RECONFIGURE_CPPUNIT)
#	set(PATCH_COMMAND_RECONFIGURE autoreconf -fvi) #autoreconf required on osx - based on contents of portfile)
#else()
#	set(PATCH_COMMAND_RECONFIGURE echo noreconf)
#endif()
#
#ExternalProject_Add(CppUnit
#	INSTALL_DIR			${HEMELB_DEPENDENCIES_INSTALL_PATH}
#	URL					${CPPUNIT_TARBALL}
#	CONFIGURE_COMMAND	<SOURCE_DIR>/configure --prefix=<INSTALL_DIR> --mandir=<INSTALL_DIR>/share/man --disable-doxygen --disable-dot "${CPPUNIT_CONFIGURE_OPTIONS}"
#	BUILD_COMMAND		make -j${HEMELB_SUBPROJECT_MAKE_JOBS}
#	PATCH_COMMAND		${PATCH_COMMAND_LDL} && ${PATCH_COMMAND_RECONFIGURE}
#	BUILD_IN_SOURCE		1)

# GOOGLE C++ TEMPLATE
# -------------------
option(CTEMPLATE_PATCH_VACOPY "Define va_copy macro through patch" OFF)
option(CTEMPLATE_PATCH_ALIGN "Define GTL align macros as gnu" OFF)

message("Building CTEMPLATE")

find_file(CTEMPLATE_TARBALL ctemplate-2.3.tar.gz
	DOC "Path to download CTEMPLATE (can be url http://)"
	PATHS ${HEMELB_DEPENDENCIES_PATH}/src)

if(NOT CTEMPLATE_TARBALL)
	message("No CTEMPLATE source found, will download")
	set(CTEMPLATE_TARBALL http://github.com/OlafvdSpek/ctemplate/archive/ctemplate-2.3.tar.gz
		CACHE STRING "Path to download CTEMPLATE (can be local file://)" FORCE)
endif()

if (CTEMPLATE_PATCH_VACOPY)
	set(PATCH_COMMAND_VACOPY patch -p1 < ${HEMELB_DEPENDENCIES_PATH}/patches/ctemplate_vacopy.diff)
else()
	set(PATCH_COMMAND_VACOPY echo novacopy)
endif()

if (CTEMPLATE_PATCH_ALIGN)
	set(PATCH_COMMAND_ALIGN patch -p1 < ${HEMELB_DEPENDENCIES_PATH}/patches/ctemplate_align.diff)
else()
	set(PATCH_COMMAND_ALIGN echo noalign)
endif()

ExternalProject_Add(CTemplate
	INSTALL_DIR			${HEMELB_DEPENDENCIES_INSTALL_PATH}
	URL					${CTEMPLATE_TARBALL}
	CONFIGURE_COMMAND	<SOURCE_DIR>/configure --prefix=<INSTALL_DIR> ${CTEMPLATE_CONFIGURE_OPTIONS}
	BUILD_COMMAND		make -j${HEMELB_SUBPROJECT_MAKE_JOBS}
	PATCH_COMMAND		${PATCH_COMMAND_ALIGN} && ${PATCH_COMMAND_VACOPY}
	BUILD_IN_SOURCE		1)

# BOOST
# -----
find_package(Boost 1.54)

if (Boost_FOUND)
	message("BOOST >= 1.54 already installed, no need to download")
else()
	message("BOOST not installed, will copy headers")

	find_file(BOOST_TARBALL boost_1_63_0.tar.gz
		DOC "Path to download BOOST (can be url http://)"
		PATHS ${HEMELB_DEPENDENCIES_PATH}/src)

	if(NOT BOOST_TARBALL)
		message("No BOOST source found, will download")
		set(BOOST_TARBALL http://sourceforge.net/projects/boost/files/boost/1.60.0/boost_1_60_0.tar.gz
			CACHE STRING "Path to download BOOST (can be local file://)" FORCE)
	endif()

	ExternalProject_Add(Boost
		INSTALL_DIR			${HEMELB_DEPENDENCIES_INSTALL_PATH}
		URL					${BOOST_TARBALL}
		CONFIGURE_COMMAND	""
		BUILD_COMMAND		""
		INSTALL_COMMAND		""
		COMMAND				${CMAKE_COMMAND} -E make_directory <INSTALL_DIR>/include/boost
		COMMAND				${CMAKE_COMMAND} -E copy_directory <SOURCE_DIR>/boost <INSTALL_DIR>/include/boost)
endif()

# ZLIB
# ----
#http://sourceforge.net/projects/libpng/files/zlib/1.2.6/zlib-1.2.6.tar.gz
find_package(ZLIB)

if (ZLIB_FOUND)
	message("ZLIB already installed, no need to download")
else()
	message("ZLIB not installed, will build from source")

	find_file(ZLIB_TARBALL zlib-1.2.6.tar.gz
		DOC "Path to download ZLIB (can be url http://)"
		PATHS ${HEMELB_DEPENDENCIES_PATH}/src)

	if(NOT ZLIB_TARBALL)
		message("No ZLIB source found, will download")
		set(ZLIB_TARBALL http://sourceforge.net/projects/libpng/files/zlib/1.2.6/zlib-1.2.6.tar.gz
			CACHE STRING "Path to download ZLIB (can be local file://)" FORCE)
	endif()

	ExternalProject_Add(zlib
		INSTALL_DIR			${HEMELB_DEPENDENCIES_INSTALL_PATH}
		URL					${ZLIB_TARBALL}
		CONFIGURE_COMMAND	<SOURCE_DIR>/configure --prefix=<INSTALL_DIR>
		BUILD_COMMAND		make -j${HEMELB_SUBPROJECT_MAKE_JOBS}
		BUILD_IN_SOURCE		1)
endif()

# TIRPC
# ----
#https://sourceforge.net/projects/libtirpc/files/libtirpc/1.2.6/libtirpc-1.2.6.tar.bz2
find_package(TIRPC REQUIRED)

if (TIRPC_FOUND)
	message("TIRPC already installed, no need to download")
else()
	message("TIRPC not installed, will build from source")

	find_file(TIRPC_TARBALL libtirpc-1.2.6.tar.bz2
		DOC "Path to download TIRPC (can be url http://)"
		PATHS ${HEMELB_DEPENDENCIES_PATH}/src)

	if(NOT TIRPC_TARBALL)
		message("No TIRPC source found, will download")
		set(TIRPC_TARBALL https://sourceforge.net/projects/libtirpc/files/libtirpc/1.2.6/libtirpc-1.2.6.tar.bz2
			CACHE STRING "Path to download TIRPC (can be local file://)" FORCE)
	endif()

	ExternalProject_Add(TIRPC
		INSTALL_DIR			${HEMELB_DEPENDENCIES_INSTALL_PATH}
		URL					${TIRPC_TARBALL}
		CONFIGURE_COMMAND	<SOURCE_DIR>/configure --prefix=<INSTALL_DIR> --disable-static --disable-gssapi
		BUILD_COMMAND		make -j${HEMELB_SUBPROJECT_MAKE_JOBS}
		BUILD_IN_SOURCE		1)
endif()

set(CMAKE_VERBOSE_MAKEFILE off)
