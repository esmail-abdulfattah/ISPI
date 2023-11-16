
#FILE(WRITE ${CMAKE_BINARY_DIR}/sympack.mak "TEST" )
#FILE(GENERATE OUTPUT ${CMAKE_BINARY_DIR}/sympack.mak INPUT ${CMAKE_BINARY_DIR}/sympack.mak_at_configured )
file(STRINGS /home/abdulfe/donttouch/pexsi_sympack_last_test/symPACK/build/sympack.mak_at_generated configfile)
string(REPLACE "-I/home/abdulfe/donttouch/pexsi_sympack_last_test/symPACK/build" "" newconfigfile "${configfile}")
string(REPLACE "/home/abdulfe/donttouch/pexsi_sympack_last_test/symPACK" "\${SYMPACK_DIR}" newconfigfile "${newconfigfile}")
string(REGEX REPLACE "[;]+" "\n" newconfigfile "${newconfigfile}")
file(WRITE /home/abdulfe/donttouch/pexsi_sympack_last_test/symPACK/build/install/include/sympack.mak "${newconfigfile}")
