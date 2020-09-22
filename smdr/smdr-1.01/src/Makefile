####################################################################
# 
#    Makefile for SMDR
#
########################## COMPILER ################################
# Uncomment the choices appropriate for your computing environment.
# (Not guaranteed to be complete. Please send info to the authors if
# you succeed in making things work on a new system.)
#
# Intel C compiler:
# =================
# CC	 	= icc
# SMDR_OPT 	= -O3 -unroll -w
#
# GNU C Compiler:
# ===============
 CC		= gcc
 SMDR_OPT 	= # -Wall
#
##### USERS SHOULD NOT NEED TO MODIFY ANYTHING BELOW THIS LINE #####

INC_DIRS = -I$(TSIL_DIR) -I$(TVIL_DIR)
LIB_DIRS = -L$(SMDR_ROOT) -L$(TSIL_DIR) -L$(TVIL_DIR)

CFLAGS	= $(INC_DIRS) $(SMDR_OPT)
LDFLAGS	= $(LIB_DIRS)
LIBS    = -lsmdr -l3vil -ltsil -lm
ARFLAGS = -rus

SMDR_SRCS = io.c \
            arguments.c \
	    messages.c \
	    utilities.c \
            fit_obs.c \
	    betas.c \
	    RGrun.c \
	    rk6_RG.c \
	    effpot.c \
	    GFermi.c \
	    Mh.c \
	    MW.c \
	    MZ.c \
	    Mt.c \
            Mlight.c \
            decouple_thZW.c \
            QCDQED_match.c \
            betas_QCDQED.c \
            RGrun_QCDQED.c

SMDR_OBJS = $(SMDR_SRCS:.c=.o)

MAKECOEFFS_OBJ = make_coeffs.o

all:	$(SMDR_ROOT)/libsmdr.a make_coeffs

$(SMDR_ROOT)/libsmdr.a: $(SMDR_OBJS)
	$(AR) $(ARFLAGS) $@ $?

make_coeffs: $(MAKECOEFFS_OBJ) $(SMDR_ROOT)/libsmdr.a
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(MAKECOEFFS_OBJ) $(LIBS)

AUTOrefmodel.h:  $(SMDR_ROOT)/benchmark_data/ReferenceModel.dat
	sed /^#/d < $< > $@

clean:
	rm -rf libsmdr.a core* a.out make_coeffs *~ *.o *.T 

# ===== Header and other dependencies: =====
# This could be more selective, with only specific objs that require
# tsil and/or 3vil given in each case:
$(SMDR_OBJS) MAKE_COEFFS_OBJ: smdr.h smdr_internal.h \
                              $(TSIL_DIR)/tsil.h \
                              $(TVIL_DIR)/3vil.h

io.o: smdr_ioparams.h

fit_obs.o: AUTOrefmodel.h
