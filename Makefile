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


#This is a list of shared directories. The paths are relative to the directories in which the codes reside
#If new shared directories are created they should be added here
export GLOBAL_MODULES = ../Finite_diff ../Force ../Input ../MG ../Misc ../US_PP ../XC
export ON_NEGF_MODULES = ../ON/ON-NEGF-share


# This shell script will remove object files in shared directories, if they were created after the code was last compiled
# This normally means that the objects were compiled for a different code, so it is the best to force their recompile
define clean-global
	if [ -e .build.log ];\
	then for module in $(GLOBAL_MODULES);\
		do for file in $$module/*.o;\
			do if [ $$file -nt .build.log ] ;\
			   then echo "Removing $$file"; rm $$file;\
		       fi;\
			done;\
	 	 done;\
	else for module in $(GLOBAL_MODULES);\
		do for file in $$module/*.o; \
			 do echo "Removing $$file";\
				if [ -e $$file ]; \
				then rm $$file;\
				else echo "$$file not exist";\
				fi;\
			 done; \
		done; \
	fi
endef

define clean-on-negf-share
	if [ -e .build.log ];\
	then for module in $(ON_NEGF_MODULES);\
		do for file in $$module/*.o;\
			do if [ $$file -nt .build.log ] ;\
			   then echo "Removing $$file"; rm $$file;\
		       fi;\
			done;\
	 	 done;\
	else for module in $(ON_NEGF_MODULES);\
		do for file in $$module/*.o; \
			 do echo "Removing $$file";\
				if [ -e $$file ]; \
				then rm $$file;\
				else echo "$$file not exist";\
				fi;\
			 done; \
		done; \
	fi
endef

all:  rmg-linux on-linux NEGF-linux

all-xt: rmg-xt  on-xt  NEGF-xt

all-aix: rmg-aix  on-aix NEGF-aix


# RMG, real-space US global grid code:
rmg-linux: RMG/Headers/make_conf.h
	@echo "#define MPI 1" > RMG/Headers/arch.h
	@echo "#define LINUX 1" >> RMG/Headers/arch.h
	cd lib/libxc/; $(MAKE) -f Make.rmg
	@cd RMG; $(clean-global)
	@echo "#define BUILD_DATE \"`date \"+%c\"`\"" >> RMG/Headers/arch.h
	cd RMG; $(MAKE) -j 8 -f Make.linux 2>&1 | tee .build.log
	

rmg-xt: RMG/Headers/make_conf.h
	@echo "#define MPI 1" > RMG/Headers/arch.h
	@echo "#define LINUX 1" >> RMG/Headers/arch.h
	cd lib/libxc/; $(MAKE) -f Make.rmg
	@cd RMG; $(clean-global)
	@echo "#define BUILD_DATE \"`date \"+%c\"`\"" >> RMG/Headers/arch.h
	cd RMG; $(MAKE) -f Make.xt 2>&1 | tee .build.log

rmg-aix: RMG/Headers/make_conf.h
	@echo "#define MPI 1" > RMG/Headers/arch.h
	@echo '#define AIX 1' >> RMG/Headers/arch.h
	cd lib/libxc/; $(MAKE) -j 8 -f Make.rmg
	@cd RMG; $(clean-global)
	@echo "#define BUILD_DATE \"`date \"+%c\"`\"" >> RMG/Headers/arch.h
	cd RMG;gmake -f Make.aix 2>&1 | tee .build.log


# Order-N targets go here
on-linux: 
	@echo "#define LINUX 1" > ON/Headers/arch.h
	@echo "#define MPI 1" >> ON/Headers/arch.h
	@echo "#define HYBRID_MODEL 0" >> ON/Headers/arch.h;
	@echo "#define THREADS_PER_NODE 1" >> ON/Headers/arch.h;
	@cd ON; $(clean-global); $(clean-on-negf-share)
	cd ON; $(MAKE) -j 8 -f Make.linux 2>&1 | tee .build.log

on-xt: 
	@echo "#define LINUX 1" > ON/Headers/arch.h
	@echo "#define MPI 1" >> ON/Headers/arch.h
	@echo "#define HYBRID_MODEL 0" >> ON/Headers/arch.h;
	@echo "#define THREADS_PER_NODE 1" >> ON/Headers/arch.h;
	@cd ON; $(clean-global); $(clean-on-negf-share)
	cd ON; $(MAKE) -f Make.xt 2>&1 | tee .build.log


on-aix: 
	@echo '#define AIX_MPI 1' > ON/Headers/arch.h
	@echo "#define PARALLEL_MESSAGE 1" >> ON/Headers/arch.h
	@echo "#define HYBRID_MODEL 0" >> ON/Headers/arch.h;
	@echo "#define THREADS_PER_NODE 1" >> ON/Headers/arch.h;
	@cd ON; $(clean-global); $(clean-on-negf-share)
	cd ON; gmake -f Make.aix 2>&1 | tee .build.log


# Targets for NEGF compilation go here
NEGF-linux: 
	@echo '#define LINUX 1' > NEGF/Headers/arch.h
	@echo "#define MPI 1" >> NEGF/Headers/arch.h
	@echo "#define HYBRID_MODEL 0" >> NEGF/Headers/arch.h;
	@echo "#define THREADS_PER_NODE 1" >> NEGF/Headers/arch.h;
	@cd NEGF; $(clean-global); $(clean-on-negf-share)
	cd NEGF; $(MAKE) -j 8 -f Make.linux 2>&1 | tee .build.log

NEGF-xt: 
	@echo '#define LINUX 1' > NEGF/Headers/arch.h
	@echo "#define MPI 1" >> NEGF/Headers/arch.h
	@echo "#define HYBRID_MODEL 0" >> NEGF/Headers/arch.h;
	@echo "#define THREADS_PER_NODE 1" >> NEGF/Headers/arch.h;
	@cd NEGF; $(clean-global); $(clean-on-negf-share)
	cd NEGF; $(MAKE) -f Make.xt 2>&1 | tee .build.log

NEGF-aix: 
	@echo '#define AIX_MPI 1' > NEGF/Headers/arch.h
	@echo "#define PARALLEL_MESSAGE 1" >> NEGF/Headers/arch.h
	@echo "#define HYBRID_MODEL 0" >> NEGF/Headers/arch.h;
	@echo "#define THREADS_PER_NODE 1" >> NEGF/Headers/arch.h;
	@cd NEGF; $(clean-global); $(clean-on-negf-share)
	cd NEGF; gmake -f Make.aix 2>&1 | tee .build.log

#Clean targets

.PHONY: clean-common
clean-common:
	cd RMG; find $(GLOBAL_MODULES)  \( -name '*.o' -o -name '*.oo' \) -exec rm {} \;

clean-rmg: clean-common
	find RMG  \( -name '*.o' -o -name '*.oo' \) -exec rm {} \;

clean-on: clean-common
	find ON  \( -name '*.o' -o -name '*.oo' \) -exec rm {} \;

clean-NEGF: clean-common
	find NEGF ON/ON-NEGF-share  \( -name '*.o' -o -name '*.oo' \) -exec rm {} \;

clean-all: clean-common clean-rmg clean-on clean-NEGF 

clean-libxc:
	cd lib/libxc; make clean 



RMG/Headers/make_conf.h:
	@echo "/****** Compile-time options for real-space code are set here  ******/" > RMG/Headers/make_conf.h;
	@echo "" >> RMG/Headers/make_conf.h;
	@echo "/* Number of processors */" >> RMG/Headers/make_conf.h;
	@echo "#define NPES 2" >> RMG/Headers/make_conf.h;
	@echo "" >> RMG/Headers/make_conf.h;
	@echo "/* 3D processor grid PE_X*PE_Y*PE_Z must equal NPES */ " >> RMG/Headers/make_conf.h;
	@echo "#define PE_X 2" >> RMG/Headers/make_conf.h;
	@echo "#define PE_Y 1" >> RMG/Headers/make_conf.h;
	@echo "#define PE_Z 1" >> RMG/Headers/make_conf.h;
	@echo "" >> RMG/Headers/make_conf.h;
	@echo "/* To enable MPI/PThreads hybrid model, experimental for now */ " >> RMG/Headers/make_conf.h;
	@echo "#define HYBRID_MODEL 0" >> RMG/Headers/make_conf.h;
	@echo "/* Number of threads in MPI/PThreads mode */ " >> RMG/Headers/make_conf.h;
	@echo "#define THREADS_PER_NODE 4" >> RMG/Headers/make_conf.h;
	@echo "" >> RMG/Headers/make_conf.h;
	@echo "/* Gamma point only set to 1, otherwise, 0 */" >> RMG/Headers/make_conf.h;
	@echo "#define GAMMA_PT 1" >> RMG/Headers/make_conf.h;
	@echo "" >> RMG/Headers/make_conf.h;
	@echo "/* Number of points in coarse grid in x,y and z directions*/" >> RMG/Headers/make_conf.h;
	@echo "#define NX_GRID 48" >> RMG/Headers/make_conf.h;
	@echo "#define NY_GRID 48" >> RMG/Headers/make_conf.h;
	@echo "#define NZ_GRID 48" >> RMG/Headers/make_conf.h;
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
	@echo "/* Papi performance monitors */" >> RMG/Headers/make_conf.h;
	@echo "#define PAPI_PERFMON 0" >> RMG/Headers/make_conf.h;
	@echo "/* Experimental fast Mehrstellen operator. Disabled by default. */" >> RMG/Headers/make_conf.h;
	@echo "#define FAST_MEHR 1" >> RMG/Headers/make_conf.h;
	@echo "";
	@echo "/* Experimental fast ortho. Disabled by default. */" >> RMG/Headers/make_conf.h;
	@echo "#define FAST_ORTHO 0" >> RMG/Headers/make_conf.h;
	@echo "";
	@echo "/* Use GPU accelerations. */" >> RMG/Headers/make_conf.h;
	@echo "#define GPU_ENABLED 0" >> RMG/Headers/make_conf.h;
	@echo "";
	@echo "ERROR: File RMG/Headers/make_conf.h does not exist"; 
	@echo "Headers/make_conf.h is set to default, configure it for your system before compiling !!! ";
	@echo "";
	exit 1; 

