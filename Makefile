#************************** SVN Revision Information **************************
#**    $Id$    **
#*****************************************************************************/
 
#
#	Top-level makefile for QMD
#
# QMD
# Version 2.0
#
# Revision 1.1  Mar, 3 2000.  Emil Briggs
#
#


SHELL = /bin/sh

all:  rmg-linux clean-common on-linux NEGF-linux

all-xt: rmg-xt clean-common on-xt clean-common NEGF-xt

all-aix: rmg-aix clean-common on-aix clean-common NEGF-aix


# RMG, real-space US global grid code:
rmg-linux: RMG/Headers/make_conf.h
	@echo "#define MPI 1" > RMG/Headers/arch.h
	@echo "#define LINUX 1" >> RMG/Headers/arch.h
	cd lib/libxc/; $(MAKE) -f Make.rmg
	cd RMG; $(MAKE) -j 8 -f Make.linux

rmg-xt: RMG/Headers/make_conf.h
	@echo "#define MPI 1" > RMG/Headers/arch.h
	@echo "#define LINUX 1" >> RMG/Headers/arch.h
	cd ../lib/libxc/; $(MAKE) -f Make.rmg
	cd RMG; $(MAKE) -f Make.xt

rmg-aix: RMG/Headers/make_conf.h
	@echo "#define MPI 1" > RMG/Headers/arch.h
	@echo '#define AIX 1' >> RMG/Headers/arch.h
	cd ../lib/libxc/; $(MAKE) -j 8 -f Make.rmg
	cd RMG;gmake -f Make.aix

# Order-N targets go here
on-linux: 
	@echo "#define LINUX 1" > ON/Headers/arch.h
	@echo "#define MPI 1" >> ON/Headers/arch.h
	cd ON; $(MAKE) -j 8 -f Make.linux

on-xt: 
	@echo "#define LINUX 1" > ON/Headers/arch.h
	@echo "#define MPI 1" >> ON/Headers/arch.h
	cd ON; $(MAKE) -f Make.xt

on-aix: 
	@echo '#define AIX_MPI 1' > ON/Headers/arch.h
	@echo "#define PARALLEL_MESSAGE 1" >> ON/Headers/arch.h
	cd ON; gmake -f Make.aix


# Targets for NEGF compilation go here
NEGF-linux: 
	@echo '#define LINUX 1' > NEGF/Headers/arch.h
	@echo "#define MPI 1" >> NEGF/Headers/arch.h
	cd NEGF; $(MAKE) -j 8 -f Make.linux

NEGF-xt: 
	@echo '#define LINUX 1' > NEGF/Headers/arch.h
	@echo "#define MPI 1" >> NEGF/Headers/arch.h
	cd NEGF; $(MAKE) -f Make.xt

NEGF-aix: 
	@echo '#define AIX_MPI 1' > NEGF/Headers/arch.h
	@echo "#define PARALLEL_MESSAGE 1" >> NEGF/Headers/arch.h
	cd NEGF; gmake -f Make.aix

#Clean targets

.PHONY: clean-common
clean-common:
	find Finite_diff/ Force/ Input/ MG/ Misc/ US_PP/ XC/  \( -name '*.o' -o -name '*.oo' \) -exec rm {} \;

clean-rmg: clean-common
	find RMG  \( -name '*.o' -o -name '*.oo' \) -exec rm {} \;

clean-on: clean-common
	find ON  \( -name '*.o' -o -name '*.oo' \) -exec rm {} \;

clean-NEGF: clean-common
	find NEGF ON/ON-NEGF-share  \( -name '*.o' -o -name '*.oo' \) -exec rm {} \;

clean-all: clean-common clean-rmg clean-on clean-NEGF 

clean-lib:
	find lib \( -name '*.o' -o -name '*.oo' \) -exec rm {} \;

clean-libxc:
	find lib/libxc \( -name '*.o' -o -name '*.oo' \) -exec rm {} \;



RMG/Headers/make_conf.h:
	@echo "/****** Compile-time options for real-space code are set here  ******/" > RMG/Headers/make_conf.h;
	@echo "" >> RMG/Headers/make_conf.h;
	@echo "/* Number of processors */" >> RMG/Headers/make_conf.h;
	@echo "#define NPES 2" >> RMG/Headers/make_conf.h;
	@echo "" >> RMG/Headers/make_conf.h;
	@echo "/* 3D processor grid PE_X*PE_Y*PE_Z must equal NPES */ " >> RMG/Headers/make_conf.h;
	@echo "#define PE_X 1" >> RMG/Headers/make_conf.h;
	@echo "#define PE_Y 1" >> RMG/Headers/make_conf.h;
	@echo "#define PE_Z 2" >> RMG/Headers/make_conf.h;
	@echo "" >> RMG/Headers/make_conf.h;
	@echo "/* Gamma point only set to 1, otherwise, 0 */" >> RMG/Headers/make_conf.h;
	@echo "#define GAMMA_PT 1" >> RMG/Headers/make_conf.h;
	@echo "" >> RMG/Headers/make_conf.h;
	@echo "/* Number of points in coarse grid in x,y and z directions*/" >> RMG/Headers/make_conf.h;
	@echo "#define NX_GRID 32" >> RMG/Headers/make_conf.h;
	@echo "#define NY_GRID 32" >> RMG/Headers/make_conf.h;
	@echo "#define NZ_GRID 32" >> RMG/Headers/make_conf.h;
	@echo "" >> RMG/Headers/make_conf.h;
	@echo "/* How many times the fine grid is finer than the coarse grid " >> RMG/Headers/make_conf.h;
	@echo " * All three numbers have to be the same */" >> RMG/Headers/make_conf.h;
	@echo "#define FG_NX 2" >> RMG/Headers/make_conf.h;
	@echo "#define FG_NY 2" >> RMG/Headers/make_conf.h;
	@echo "#define FG_NZ 2" >> RMG/Headers/make_conf.h;
	@echo "" >> RMG/Headers/make_conf.h;
	@echo "/* Set this to 0 to turn off memory Smart-ALLOCation. (salloc.c, salloc.h)" >> RMG/Headers/make_conf.h;
	@echo " * (there is no significant time improvement)*/" >> RMG/Headers/make_conf.h;
	@echo "#define USE_SALLOC 1" >> RMG/Headers/make_conf.h;
	@echo "" >> RMG/Headers/make_conf.h;
	@echo "/* Set this to 1 if you want to use finite difference method for calculating" >> RMG/Headers/make_conf.h;
	@echo " * derivatives of beta. This is faster since it avoids doing 3 backwards fourier" >> RMG/Headers/make_conf.h;
	@echo " * transforms per ion, but it may not be very accurate since the finite diff" >> RMG/Headers/make_conf.h;
	@echo " * derivative is done on the coarse grid." >> RMG/Headers/make_conf.h;
	@echo " * Leave this set to 0 unless you know what you are doing */" >> RMG/Headers/make_conf.h;
	@echo "#define FDIFF_BETA 0" >> RMG/Headers/make_conf.h;
	@echo "" >> RMG/Headers/make_conf.h;
	@echo "/* Extra fine timers, may cause some slowdown, but they are generally useful */" >> RMG/Headers/make_conf.h;
	@echo "#define MD_TIMERS 1" >> RMG/Headers/make_conf.h;
	@echo "";
	@echo "ERROR: File Headers/make_conf.h does not exist"; 
	@echo "Headers/make_conf.h is set to default, configure it for your system before compiling !!! ";
	@echo "";
	exit 1; 

