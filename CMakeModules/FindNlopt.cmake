################################################################################
# Module to find libnlopt                                                      #
#                                                                              #
# This sets the following variables:                                           #
#   - Nlopt_FOUND                                                              #
#   - NLOPT_LIBRARIES                                                          #
#   - NLOPT_INCLUDE_PATH                                                       #
#                                                                              #
################################################################################

IF (NOT Nlopt_FOUND)

FIND_LIBRARY (NLOPT_LIBRARIES nlopt HINTS $ENV{NLOPTROOT}/lib)

FIND_PATH (NLOPT_INCLUDE_PATH nlopt.h HINTS $ENV{NLOPTROOT}/include)

ENDIF (NOT Nlopt_FOUND)
