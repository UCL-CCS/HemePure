# check for header file
find_path(CPPUNIT_INCLUDE_DIRS cppunit/Test.h
	HINTS ${CPPUNIT_DIR}/include $ENV{CPPUNIT_DIR}/include
	PATH_SUFFIXES cppunit
	DOC "directory where the CPPUNIT header is located")

# check for cppunit
find_library(CPPUNIT_LIBRARY
	NAMES cppunit
	HINTS ${CPPUNIT_DIR}/lib $ENV{CPPUNIT_DIR}/lib
	NO_DEFAULT_PATH
	DOC "the CPPUNIT library")
find_library(CPPUNIT_LIBRARY
	NAMES cppunit
	DOC "the CPPUNIT library")

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(CPPUNIT
	CPPUNIT_INCLUDE_DIR CPPUNIT_LIBRARY)
