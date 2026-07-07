#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "DG_ExProj::DG_ExProj" for configuration "Release"
set_property(TARGET DG_ExProj::DG_ExProj APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(DG_ExProj::DG_ExProj PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libdg_exproj.so"
  IMPORTED_SONAME_RELEASE "libdg_exproj.so"
  )

list(APPEND _IMPORT_CHECK_TARGETS DG_ExProj::DG_ExProj )
list(APPEND _IMPORT_CHECK_FILES_FOR_DG_ExProj::DG_ExProj "${_IMPORT_PREFIX}/lib/libdg_exproj.so" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
