cmake_path(GET CMAKE_CURRENT_SOURCE_DIR STEM DIRNAME)
set(PROJECT_PREFIX ${PROJECT_PREFIX}${DIRNAME})

add_subdirectory(AnisotropicCompartments)
add_subdirectory(ODF)
add_subdirectory(UnweightedMCM)
add_subdirectory(WeightedMCM)
