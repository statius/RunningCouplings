####################################################################
# 
#    Makefile for SMDR APPLICATIONS
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

INC_DIRS = -I$(SMDR_ROOT)
LIB_DIRS = -L$(SMDR_ROOT)

CFLAGS	= $(INC_DIRS) $(SMDR_OPT)
LDFLAGS	= $(LIB_DIRS)
LIBS    = -lsmdr -l3vil -ltsil -lm

APP_SRCS = calc_RGrun.c \
           calc_all.c \
           calc_fit.c \
           calc_m2.c \
           calc_vev.c \
           fig_GFermi_vs_Q.c \
           fig_MW_vs_Q.c \
           fig_MZ_vs_Q.c \
           fig_Mh_vs_Q.c \
           fig_Mt_vs_Q.c \
           fig_RGrun_QCDQED.c \
           fig_RGrun_vs_Q.c \
           fig_lambda_m2_vs_Q.c \
           fig_lambda_vs_Q.c \
           fig_m2_vs_Q.c \
           fig_vev_vs_Q.c

# This works fine too!
# APP_SRCS = $(shell ls *.c)

APP_NAMES = $(APP_SRCS:.c=)

APP_OBJS = $(APP_SRCS:.c=.o)

SMDR_APPS = $(APP_NAMES:%=$(SMDR_ROOT)/%)

all:    $(SMDR_APPS) benchmarkdata

$(SMDR_APPS): $(SMDR_ROOT)/%: %.o \
              $(SMDR_ROOT)/libsmdr.a \
              $(SMDR_ROOT)/libtsil.a \
              $(SMDR_ROOT)/lib3vil.a
	$(CC) $(CFLAGS) $(LDFLAGS) $< -o $@ $(LIBS)

benchmarkdata:
	diff $(SMDR_ROOT)/ReferenceModel.dat $(SMDR_ROOT)/benchmark_data/ReferenceModel.dat > /dev/null || cp $(SMDR_ROOT)/benchmark_data/ReferenceModel.dat $(SMDR_ROOT)

clean:
	rm -rf core* a.out $(SMDR_APPS) *~ *.o *.T

# ===== Header dependencies: =====
$(APP_OBJS):  $(SMDR_ROOT)/smdr.h \
              $(SMDR_ROOT)/tsil.h \
              $(SMDR_ROOT)/3vil.h
