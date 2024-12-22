# - Locate DMG4 library
# Defines:
#
#  DMG4_FOUND
#  DMG4L_INCLUDE_DIR
#  DMG4_DARKMATTER_INCLUDE_DIR
#  DMG4_INCLUDE_DIRS (not cached)
#  DMG4_LIBRARY
#  DMG4_DARKMATTER_LIBRARY
#  DMG4_UTILSDM_LIBRARY
#  DMG4_LIBRARIES (not cached)

find_path(DMG4_INCLUDE_DIR DarkMatterPhysics.hh
          HINTS $ENV{DMG4_ROOT_DIR}/include/DMG4 ${DMG4_ROOT_DIR}/include/DMG4)
find_path(DMG4_DARKMATTER_INCLUDE_DIR DarkMatter.hh
          HINTS $ENV{DMG4_ROOT_DIR}/include ${DMG4_ROOT_DIR}/include)
find_path(DMG4_UTILSDM_INCLUDE_DIR DarkMatterParametersRegistry.hh
          HINTS $ENV{DMG4_ROOT_DIR}/include/UtilsDM ${DMG4_ROOT_DIR}/include/UtilsDM)
#find_path(DMG4_DMPARTICLES_INCLUDE_DIR DMParticleAPrime.hh
#          HINTS $ENV{DMG4_ROOT_DIR}/include/DMG4/DMParticles ${DMG4_ROOT_DIR}/include/DMG4/DMParticles)
#find_path(DMG4_DMPROCESSES_INCLUDE_DIR DMProcessDMBrem.hh
#          HINTS $ENV{DMG4_ROOT_DIR}/include/DMG4/DMProcesses ${DMG4_ROOT_DIR}/include/DMG4/DMProcesses)

find_library(DMG4_LIBRARY NAMES DMG4
             HINTS $ENV{DMG4_ROOT_DIR}/lib ${DMG4_ROOT_DIR}/lib)
find_library(DMG4_DARKMATTER_LIBRARY NAMES DarkMatter
             HINTS $ENV{DMG4_ROOT_DIR}/lib ${DMG4_ROOT_DIR}/lib)
find_library(DMG4_UTILSDM_LIBRARY NAMES UtilsDM
             HINTS $ENV{DMG4_ROOT_DIR}/lib ${DMG4_ROOT_DIR}/lib)

set(DMG4_LIBRARIES ${DMG4_LIBRARY} ${DMG4_DARKMATTER_LIBRARY} ${DMG4_UTILSDM_LIBRARY})

#set(DMG4_INCLUDE_DIRS ${DMG4_INCLUDE_DIR} ${DMG4_DMPARTICLES_INCLUDE_DIR} ${DMG4_DMPROCESSES_INCLUDE_DIR} ${DMG4_DARKMATTER_INCLUDE_DIR})
set(DMG4_INCLUDE_DIRS ${DMG4_INCLUDE_DIR} ${DMG4_INCLUDE_DIR}/DMParticles ${DMG4_INCLUDE_DIR}/DMProcesses ${DMG4_DARKMATTER_INCLUDE_DIR} ${DMG4_UTILSDM_INCLUDE_DIR})

# handle the QUIETLY and REQUIRED arguments and set DMG4_FOUND to TRUE if
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(DMG4 DEFAULT_MSG DMG4_INCLUDE_DIRS DMG4_LIBRARIES)

mark_as_advanced(DMG4_FOUND DMG4_INCLUDE_DIR DMG4_DARKMATTER_INCLUDE_DIR DMG4_LIBRARY DMG4_DARKMATTER_LIBRARY DMG4_UTILSDM_LIBRARY)
