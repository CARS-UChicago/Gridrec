TOP=../..
include $(TOP)/configure/CONFIG
#----------------------------------------
#  ADD MACRO DEFINITIONS AFTER THIS LINE
#=============================

SHARED_LIBRARIES = YES
STATIC_BUILD = YES
LIBRARY_HOST += GridrecIDL
USR_CFLAGS += -DIDL
GridrecIDL_LDFLAGS_Linux += $(STATIC_LDFLAGS)

GridrecIDL_SRCS += grid.c GridrecIDL.c grid_math.c pswf.c filters.c fft_fftw.c

GridrecIDL_LIBS_WIN32 += libfftw3f-3
# To use the version of fftw3f included with tomoRecon use this line
GridrecIDL_LIBS_Linux += fftw3f
# To use the system version of fftw3f use this line
#tomoRecon_SYS_LIBS_Linux += fftw3f
GridrecIDL_SYS_LIBS_Darwin += fftw3f

#=============================


#===========================

include $(TOP)/configure/RULES
#----------------------------------------
#  ADD RULES AFTER THIS LINE
