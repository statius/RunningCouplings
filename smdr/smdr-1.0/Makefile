####################################################################
# 
#    TOP LEVEL Makefile for SMDR
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
################ INSTALLATION DIRECTORIES ##########################
# Set these to the locations where you would like to install the
# header files and static libraries, respectively.
#
INSTALL_INCS = /usr/include
INSTALL_LIBS = /usr/lib
#
#
##### USERS SHOULD NOT NEED TO MODIFY ANYTHING BELOW THIS LINE #####

export SMDR_ROOT = $(shell pwd)
export TSIL_DIR = $(SMDR_ROOT)/$(shell ls tsil*.tar.gz | sed -e "s/.tar.gz//")
export TVIL_DIR = $(SMDR_ROOT)/$(shell ls 3vil*.tar.gz | sed -e "s/.tar.gz//")
export SMDR_SRC = $(SMDR_ROOT)/src
export SMDR_APP_SRC = $(SMDR_ROOT)/applications

all:	libs apps

libs:	tsil 3vil smdr 

tsil:
	test -s $(TSIL_DIR)/tsil.h || gunzip -c tsil*.tar.gz | tar xvf - ;
	+$(MAKE) -C $(TSIL_DIR)

3vil:
	test -s $(TVIL_DIR)/3vil.h || gunzip -c 3vil*.tar.gz | tar xvf - ;
	+$(MAKE) -C $(TVIL_DIR)

smdr:   $(SMDR_ROOT)/libtsil.a $(SMDR_ROOT)/lib3vil.a \
        $(SMDR_ROOT)/tsil.h $(SMDR_ROOT)/3vil.h $(SMDR_ROOT)/smdr.h 
	+$(MAKE) -C $(SMDR_SRC)

apps:	$(SMDR_ROOT)/libtsil.a $(SMDR_ROOT)/lib3vil.a $(SMDR_ROOT)/libsmdr.a \
        $(SMDR_ROOT)/tsil.h $(SMDR_ROOT)/3vil.h $(SMDR_ROOT)/smdr.h 
	+$(MAKE) -C $(SMDR_APP_SRC)

$(SMDR_ROOT)/libtsil.a: $(TSIL_DIR)/libtsil.a
	cp $(TSIL_DIR)/libtsil.a $(SMDR_ROOT)/libtsil.a

$(SMDR_ROOT)/tsil.h: $(TSIL_DIR)/tsil.h
	cp $(TSIL_DIR)/tsil.h $(SMDR_ROOT)/tsil.h

$(SMDR_ROOT)/lib3vil.a: $(TVIL_DIR)/lib3vil.a
	cp $(TVIL_DIR)/lib3vil.a $(SMDR_ROOT)/lib3vil.a

$(SMDR_ROOT)/3vil.h: $(TVIL_DIR)/3vil.h
	cp $(TVIL_DIR)/3vil.h $(SMDR_ROOT)/3vil.h

$(SMDR_ROOT)/smdr.h: $(SMDR_SRC)/smdr.h
	cp $(SMDR_SRC)/smdr.h $(SMDR_ROOT)

clean:
	+$(MAKE) cleantsil
	+$(MAKE) clean3vil
	+$(MAKE) cleansmdr

cleantsil:
	+$(MAKE) clean -C $(TSIL_DIR)
	rm $(SMDR_ROOT)/libtsil.a $(SMDR_ROOT)/tsil.h

clean3vil:
	+$(MAKE) clean -C $(TVIL_DIR)
	rm $(SMDR_ROOT)/lib3vil.a $(SMDR_ROOT)/3vil.h

cleansmdr:
	+$(MAKE) clean -C $(SMDR_SRC)
	rm $(SMDR_ROOT)/libsmdr.a $(SMDR_ROOT)/smdr.h
	$(MAKE) cleanapps

cleanapps:
	+$(MAKE) clean -C $(SMDR_APP_SRC)
	rm $(SMDR_ROOT)/ReferenceModel.dat

tidy:
	+$(MAKE) clean -C $(TSIL_DIR)
	+$(MAKE) clean -C $(TVIL_DIR)
	+$(MAKE) clean -C $(SMDR_SRC)
	+$(MAKE) clean -C $(SMDR_APP_SRC)

install:   $(SMDR_ROOT)/libsmdr.a \
           $(SMDR_ROOT)/libtsil.a \
           $(SMDR_ROOT)/lib3vil.a \
           $(SMDR_ROOT)/smdr.h \
           $(SMDR_ROOT)/tsil.h \
           $(SMDR_ROOT)/3vil.h
	cp $(SMDR_ROOT)/lib*.a $(INSTALL_LIBS)
	cp $(SMDR_ROOT)/*.h $(INSTALL_INCS)
