#----------------------------------------------------------------
# Generated CMake target import file for configuration "opt".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "symPACK::symPACK" for configuration "opt"
set_property(TARGET symPACK::symPACK APPEND PROPERTY IMPORTED_CONFIGURATIONS OPT)
set_target_properties(symPACK::symPACK PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_OPT "CXX;Fortran"
  IMPORTED_LOCATION_OPT "${_IMPORT_PREFIX}/lib/libsympack.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS symPACK::symPACK )
list(APPEND _IMPORT_CHECK_FILES_FOR_symPACK::symPACK "${_IMPORT_PREFIX}/lib/libsympack.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
