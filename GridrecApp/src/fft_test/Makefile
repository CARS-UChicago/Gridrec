TOP=../../..
include $(TOP)/configure/CONFIG
#----------------------------------------
#  ADD MACRO DEFINITIONS AFTER THIS LINE
#=============================

SHARED_LIBRARIES = YES
STATIC_BUILD = YES
LIBRARY_HOST += fftTestIDL
fftTestIDL_LDFLAGS_Linux += $(STATIC_LDFLAGS)

fftTestIDL_SRCS += fft_test.c fft.c

fftTestIDL_LIBS_WIN32 += libfftw3f-3
# To use the version of fftw3f included with tomoRecon use this line
fftTestIDL_LIBS_Linux += fftw3f
# To use the system version of fftw3f use this line
#tomoRecon_SYS_LIBS_Linux += fftw3f
fftTestIDL_SYS_LIBS_Darwin += fftw3f

#=============================


#===========================

include $(TOP)/configure/RULES
#----------------------------------------
#  ADD RULES AFTER THIS LINE
