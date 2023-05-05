include(CheckCXXSourceCompiles)
include(CheckCXXSymbolExists)

CHECK_CXX_SOURCE_COMPILES("#include <sys/time.h>\n#include <sys/resource.h>\nint main(int c,char** v){ rusage usage;\ngetrusage(RUSAGE_SELF, &usage);\nreturn usage.ru_maxrss; }" HAVE_RUSAGE)

# cstdint is the c++11 version of C99 stdint.h.
# better to go with cstdint, but stdint.h is available more widely.
find_path(HAVE_STDINT_H stdint.h)
find_path(HAVE_CSTDINT cstdint)
if(HAVE_CSTDINT)
	add_definitions(-DHEMELB_HAVE_CSTDINT=TRUE)
elseif(NOT HAVE_STDINT_H)
	message(ERROR "Neither cstdint nor stdint.h found")
endif()
CHECK_CXX_SYMBOL_EXISTS("funopen" "stdio.h" HAVE_FUNOPEN)
if(HAVE_FUNOPEN)
  add_definitions(-DHAVE_FUNOPEN)
endif()
CHECK_CXX_SYMBOL_EXISTS("fopencookie" "stdio.h" HAVE_FOPENCOOKIE)
if(HAVE_FOPENCOOKIE)
  add_definitions(-DHAVE_FOPENCOOKIE)
endif()
