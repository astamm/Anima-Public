cmake_path(GET CMAKE_CURRENT_SOURCE_DIR STEM DIRNAME)
set(PROJECT_PREFIX ${PROJECT_PREFIX}${DIRNAME})

add_subdirectory(MCM)
add_subdirectory(NoiseVariance)
add_subdirectory(ODF)
add_subdirectory(Tensors)
add_subdirectory(TOD)
