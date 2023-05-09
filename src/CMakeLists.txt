cmake_minimum_required(VERSION 3.1)

project(HemeLB)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
#set(CMAKE_CXX_EXTENSIONS OFF)

# conditional compile options
# ---------------------------
option(HEMELB_BUILD_TESTS_UNIT "build the unit tests" OFF)
option(HEMELB_DEPENDENCIES_SET_RPATH "set runtime RPATH" ON)
option(HEMELB_IMAGES_TO_NULL "write images to null" ON)
option(HEMELB_STATIC_ASSERT "use simple compile-time assertions" ON)
option(HEMELB_TRACER_PARTICLES "use particles as tracers" ON)
option(HEMELB_USE_ALL_WARNINGS_GNU "show all compiler warnings on development builds (gnu-style-compilers)" OFF)
option(HEMELB_USE_SSE3 "use SSE3 intrinsics" ON)
option(HEMELB_USE_AVX2 "use AVX2 intrinsics" OFF)
option(HEMELB_USE_AVX512 "use AVX512 intrinsics" OFF)
option(HEMELB_USE_VELOCITY_WEIGHTS_FILE "use velocity weights file" OFF)
option(HEMELB_USE_GMYPLUS "use GMY+ format" OFF)
option(HEMELB_USE_MPI_CALL "use standard MPI functions when reading blocks" OFF)
option(HEMELB_USE_MPI_WIN "use MPI windowing to help load large domains" OFF)
option(HEMELB_USE_BIGMPI "use BigMPI to read data from file" OFF)
option(HEMELB_USE_PARMETIS "use ParMETIS" OFF)
#option(HEMELB_USE_ZOLTAN "use Zoltan" OFF)
#option(HEMELB_VALIDATE_GEOMETRY "validate geometry" OFF)

set(HEMELB_EXECUTABLE "hemepure"
	CACHE STRING "name of executable to produce")
set(HEMELB_READING_GROUP_SIZE 2
	CACHE INTEGER "number of cores to use to read geometry file")
set(HEMELB_READING_GROUP_SPACING 1
	CACHE INTEGER "spacing between cores to use to read geometry file")
set(HEMELB_LOG_LEVEL Info
	CACHE STRING "log level: 'Critical','Error','Warning','Info','Debug','Trace'")
set(HEMELB_DEPENDENCIES_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../dep/"
	CACHE FILEPATH "path to find dependency find modules")
set(HEMELB_DEPENDENCIES_INSTALL_PATH "${HEMELB_DEPENDENCIES_PATH}/install"
	CACHE FILEPATH "path to find dependency includes and libraries")
set(HEMELB_OPTIMISATION "-O3"
	CACHE STRING "optimisation level: <blank> or -O1 to -O4)")
set(HEMELB_LATTICE "D3Q19"
	CACHE STRING "lattice type: 'D3Q15','D3Q19','D3Q27','D3Q15i'")
set(HEMELB_KERNEL "LBGK"
	CACHE STRING "kernel: 'LBGK','EntropicAnsumali','EntropicChik','MRT','TRT','NNCY','NNCYMOUSE','NNC','NNTPL'")
set(HEMELB_WALL_BOUNDARY "BFL"
	CACHE STRING "boundary conditions to be used at the walls: 'BFL','GZS','SIMPLEBOUNCEBACK','JUNKYANG'")
set(HEMELB_INLET_BOUNDARY "NASHZEROTHORDERPRESSUREIOLET"
	CACHE STRING "boundary conditions to be used at the inlet: 'NASHZEROTHORDERPRESSUREIOLET','LADDIOLET'")
set(HEMELB_OUTLET_BOUNDARY "NASHZEROTHORDERPRESSUREIOLET"
	CACHE STRING "boundary conditions to be used at the outlets: 'NASHZEROTHORDERPRESSUREIOLET','LADDIOLET'")
set(HEMELB_COMPUTE_ARCHITECTURE "NEUTRAL"
	CACHE STRING "architecture of the machine being used 'INTELSANDYBRIDGE','AMDBULLDOZER','NEUTRAL','ISBFILEVELOCITYINLET'")
set(HEMELB_WALL_INLET_BOUNDARY "NASHZEROTHORDERPRESSUREBFL"
	CACHE STRING "boundary conditions to be used at corners between walls and inlets: 'NASHZEROTHORDERPRESSURESBB','NASHZEROTHORDERPRESSUREBFL','LADDIOLETSBB','LADDIOLETBFL'")
set(HEMELB_WALL_OUTLET_BOUNDARY "NASHZEROTHORDERPRESSUREBFL"
	CACHE STRING "boundary conditions to be used at corners between walls and outlets: 'NASHZEROTHORDERPRESSURESBB','NASHZEROTHORDERPRESSUREBFL','LADDIOLETSBB','LADDIOLETBFL'")
set(HEMELB_POINTPOINT_IMPLEMENTATION Coalesce
	CACHE STRING "point to point comms implementation: 'Coalesce','Separated','Immediate'" )
set(HEMELB_GATHERS_IMPLEMENTATION Separated
	CACHE STRING "gather comms implementation: 'Separated','ViaPointPoint'" )
set(HEMELB_ALLTOALL_IMPLEMENTATION Separated
	CACHE STRING "alltoall comms implementation: 'Separated','ViaPointPoint'" )
option(HEMELB_SEPARATE_CONCERNS "communicate for each concern separately" OFF)

# stop warnings on ARCHER
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated-declarations")

# add warnings flags to development build types
if(HEMELB_USE_ALL_WARNINGS_GNU)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
	set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall")
endif()

if(HEMELB_OPTIMISATION)
	set(CMAKE_CXX_FLAGS_RELEASE ${HEMELB_OPTIMISATION})
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${HEMELB_OPTIMISATION}")
endif()

add_definitions(-DHEMELB_CODE)
add_definitions(-DHEMELB_READING_GROUP_SIZE=${HEMELB_READING_GROUP_SIZE})
add_definitions(-DHEMELB_READING_GROUP_SPACING=${HEMELB_READING_GROUP_SPACING})
add_definitions(-DHEMELB_LATTICE=${HEMELB_LATTICE})
add_definitions(-DHEMELB_KERNEL=${HEMELB_KERNEL})
add_definitions(-DHEMELB_WALL_BOUNDARY=${HEMELB_WALL_BOUNDARY})
add_definitions(-DHEMELB_INLET_BOUNDARY=${HEMELB_INLET_BOUNDARY})
add_definitions(-DHEMELB_OUTLET_BOUNDARY=${HEMELB_OUTLET_BOUNDARY})
add_definitions(-DHEMELB_WALL_INLET_BOUNDARY=${HEMELB_WALL_INLET_BOUNDARY})
add_definitions(-DHEMELB_WALL_OUTLET_BOUNDARY=${HEMELB_WALL_OUTLET_BOUNDARY})
add_definitions(-DHEMELB_COMPUTE_ARCHITECTURE=${HEMELB_COMPUTE_ARCHITECTURE})
add_definitions(-DHEMELB_LOG_LEVEL=${HEMELB_LOG_LEVEL})

if(HEMELB_IMAGES_TO_NULL)
	add_definitions(-DHEMELB_IMAGES_TO_NULL)
endif()

if(NOT HEMELB_STATIC_ASSERT)
	add_definitions(-DHEMELB_NO_STATIC_ASSERT)
endif()

if(HEMELB_TRACER_PARTICLES)
	add_definitions(-DHEMELB_TRACER_PARTICLES)
endif()

if(HEMELB_USE_SSE3)
	add_definitions(-DHEMELB_USE_SSE3)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse3")
	set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -msse3")
	set(CMAKE_CXX_FLAGS_RELEASE "${HEMELB_OPTIMISATION} -msse3")
endif()

if(HEMELB_USE_AVX2)
        add_definitions(-DHEMELB_USE_AVX2)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=core-avx2")
        set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -march=core-avx2")
        set(CMAKE_CXX_FLAGS_RELEASE "${HEMELB_OPTIMISATION} -march=core-avx2")
endif()

if(HEMELB_USE_AVX512)
        add_definitions(-DHEMELB_USE_AVX512)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -xCORE-AVX512")
        set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -xCORE-AVX512")
        set(CMAKE_CXX_FLAGS_RELEASE "${HEMELB_OPTIMISATION} -xCORE-AVX512") #-march=skylake-avx512
endif()

if(HEMELB_USE_VELOCITY_WEIGHTS_FILE)
	add_definitions(-DHEMELB_USE_VELOCITY_WEIGHTS_FILE)
endif()

if(HEMELB_USE_PARMETIS)
	add_definitions(-DHEMELB_USE_PARMETIS)
endif()

if(HEMELB_USE_GMYPLUS)
	add_definitions(-DHEMELB_USE_GMYPLUS)
endif()

if(HEMELB_USE_MPI_CALL)
	add_definitions(-DHEMELB_USE_MPI_CALL)
endif()

if(HEMELB_USE_MPI_WIN)
	add_definitions(-DHEMELB_USE_MPI_WIN)
endif()

if(HEMELB_USE_BIGMPI)
	add_definitions(-DHEMELB_USE_BIGMPI)
endif()

#if(HEMELB_VALIDATE_GEOMETRY)
#	add_definitions(-DHEMELB_VALIDATE_GEOMETRY)
#endif()

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" "${HEMELB_DEPENDENCIES_PATH}/modules/")
list(APPEND CMAKE_INCLUDE_PATH ${HEMELB_DEPENDENCIES_INSTALL_PATH}/include)
list(APPEND CMAKE_LIBRARY_PATH ${HEMELB_DEPENDENCIES_INSTALL_PATH}/lib)

# set up runtime search path for dlls
#------------------------------------
if(HEMELB_DEPENDENCIES_SET_RPATH)
	set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
	set(CMAKE_INSTALL_RPATH "${HEMELB_DEPENDENCIES_INSTALL_PATH}/lib")
endif()

include(platform_checks)

add_definitions(-DHAVE_STD_ISNAN)	# add_definitions(-DHAVE_ISNAN) is the alternative
add_definitions(-DLINUX_SCANDIR)	# definitely using Linux!!

if(HAVE_STD_ISNAN)
	add_definitions(-DHAVE_STD_ISNAN)
endif()

if(HAVE_RUSAGE)
	add_definitions(-DHAVE_RUSAGE)
endif()

include(mpi)
include(dependencies)

#if(HEMELB_USE_ZOLTAN)
#
#	set(HEMELB_DECOMPOSITION "PARMETIS"
#		CACHE STRING "decomposition method: 'PARMETIS','PTSCOTCH','BLOCK'")
#
#	if(
#		HEMELB_DECOMPOSITION STREQUAL "PARMETIS")
#		add_definitions(-DHEMELB_USE_PARMETIS)
#	elseif(
#		HEMELB_DECOMPOSITION STREQUAL "PTSCOTCH")
#		add_definitions(-DHEMELB_USE_PTSCOTCH)
#	elseif(
#		HEMELB_DECOMPOSITION STREQUAL "BLOCK")
#		add_definitions(-DHEMELB_USE_BLOCK)
#	endif()
#
#	set(DECOMPOSE_LIBRARIES ${ZOLTAN_LIBRARY} ${PTSCOTCH_LIBRARIES} ${PARMETIS_LIBRARIES})
#	add_definitions(-DHEMELB_USE_ZOLTAN)
#
#elseif(NOT HEMELB_USE_ZOLTAN)
set(DECOMPOSE_LIBRARIES ${PARMETIS_LIBRARIES})
#endif()

# resources
# ---------
set(BUILD_RESOURCE_PATH ${PROJECT_BINARY_DIR}/resources)
file(MAKE_DIRECTORY ${BUILD_RESOURCE_PATH})
set(INSTALL_RESOURCE_PATH ${CMAKE_INSTALL_PREFIX}/share/hemelb/resources)

configure_file(
	"${PROJECT_SOURCE_DIR}/resources/path_parameters.h.in"
	"${PROJECT_BINARY_DIR}/resources/path_parameters.h")

# add the binary tree to the search path for include files
# so that we will find TutorialConfig.h
include_directories("${PROJECT_BINARY_DIR}")
include_directories("${HEMELB_DEPENDENCIES_INSTALL_PATH}/include")

# HemeLB
# ------
set(root_sources SimulationMaster.cc)
add_executable(${HEMELB_EXECUTABLE} main.cc ${root_sources})

include_directories(${PROJECT_SOURCE_DIR})
set(package_subdirs
	configuration
	extraction
	reporting
	geometry
	lb
	net
	util
	io
	log
	colloids)

foreach(subdir ${package_subdirs})
	set(lib "hemelb_${subdir}")
	list(APPEND HEMELB_LIBRARIES ${lib})
	add_subdirectory(${subdir})
endforeach()

add_subdirectory(resources)
target_link_libraries(
	${HEMELB_EXECUTABLE}
	${HEMELB_LIBRARIES}
	${MPI_LIBRARIES}
	${DECOMPOSE_LIBRARIES}
	${TINYXML_LIBRARY}
	${Boost_LIBRARIES}
	${CTEMPLATE_LIBRARY}
	${ZLIB_LIBRARIES}
	${TIRPC_LIBRARY})
install(TARGETS ${HEMELB_EXECUTABLE} RUNTIME DESTINATION bin)
list(APPEND RESOURCES resources/report.txt.ctp resources/report.xml.ctp)

# HemeLB unittests
# ----------------
if(HEMELB_BUILD_TESTS_UNIT)
	find_package(CPPUNIT REQUIRED)
	include_directories(${CPPUNIT_INCLUDE_DIR})
	add_executable(unittests_hemelb ${root_sources})
	add_subdirectory(unittests)
	target_link_libraries(unittests_hemelb
		hemelb_unittests
		${HEMELB_LIBRARIES}
		${MPI_LIBRARIES}
		${DECOMPOSE_LIBRARIES}
		${TINYXML_LIBRARY}
		${CPPUNIT_LIBRARY}
		${Boost_LIBRARIES}
		${CTEMPLATE_LIBRARY}
		${ZLIB_LIBRARIES}
		${TIRPC_LIBRARY}
		${CMAKE_DL_LIBS})
	install(TARGETS unittests_hemelb RUNTIME DESTINATION bin)
	list(APPEND RESOURCES
		unittests/resources/four_cube.gmy
		unittests/resources/four_cube.xml
		unittests/resources/config.xml unittests/resources/config0_2_0.xml
		unittests/resources/config_file_inlet.xml unittests/resources/iolet.txt
		unittests/resources/config-velocity-iolet.xml unittests/resources/config_new_velocity_inlets.xml
		unittests/resources/velocity_inlet.txt.weights.txt
		unittests/resources/xmltest.xml unittests/resources/config_file_velocity_inlet.xml unittests/resources/velocity_inlet.txt)
endif()

# copy and install resources
# --------------------------
foreach(resource ${RESOURCES})
	configure_file(${PROJECT_SOURCE_DIR}/${resource} ${BUILD_RESOURCE_PATH} COPYONLY)
	get_filename_component(resource_leaf ${resource} NAME)
	install(FILES ${BUILD_RESOURCE_PATH}/${resource_leaf} DESTINATION share/hemelb/resources)
endforeach()
