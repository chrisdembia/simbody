# Build and install TaskSpaceControl-UR10 example, from within a
# Simbody source build. This must be invoked as a CMake subdirectory from
# the main examples CMakeLists.txt file.

set(ex_name TaskSpaceControl-UR10)
set(header_files UR10.h)
set(source_files TaskSpaceControl-UR10.cpp UR10.cpp)

FILE(GLOB geodir_geom "geometry/*")

SIMBODY_EXAMPLE("${ex_name}" "${header_files}" "${source_files}"
    "${geodir_geom}")
