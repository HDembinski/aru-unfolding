################################################################################
# Module to find Root                                                          #
#                                                                              #
# This sets the following variables:                                           #
#   - Root_FOUND                                                               #
#   - ROOT_BIN_DIR                                                             #
#   - ROOT_LIBRARIES                                                           #
#   - ROOT_INCLUDE_PATH                                                        #
#                                                                              #
################################################################################

IF (NOT Root_FOUND)

  FIND_PATH (ROOT_BIN_DIR root-config HINTS $ENV{ROOTSYS}/bin)

  GET_FILENAME_COMPONENT (ROOTSYS ${ROOT_BIN_DIR} PATH)

  IF (NOT ENV{ROOTSYS})
    SET (ENV{ROOTSYS} ${ROOTSYS})
  ENDIF (NOT ENV{ROOTSYS})

  IF (ROOTSYS)
    # --------------------------------------------------------------------------
    # Set ROOT compilation flags.
    # --------------------------------------------------------------------------
    SET (Root_FOUND TRUE)

    EXECUTE_PROCESS (COMMAND ${ROOT_BIN_DIR}/root-config --libs
      OUTPUT_VARIABLE ROOT_LIBRARIES
      OUTPUT_STRIP_TRAILING_WHITESPACE
    )

    EXECUTE_PROCESS (COMMAND ${ROOT_BIN_DIR}/root-config --incdir
      OUTPUT_VARIABLE ROOT_INCLUDE_PATH
      OUTPUT_STRIP_TRAILING_WHITESPACE
    )
  ENDIF (ROOTSYS)

ENDIF (NOT Root_FOUND)
