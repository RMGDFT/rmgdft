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

codes_dir = ~/SVN/codes
RMG_dir = $(codes_dir)/RMG
NEGF_dir = $(codes_dir)/NEGF
ON_dir = $(codes_dir)/TDDFT

#This is a list of shared directories. The paths are relative to the directories in which the codes reside
#If new shared directories are created they should be added here
export GLOBAL_MODULES = $(codes_dir)/Finite_diff $(codes_dir)/Force $(codes_dir)/Input $(codes_dir)/MG $(codes_dir)/Misc $(codes_dir)/US_PP $(codes_dir)/XC
export ON_NEGF_MODULES = $(codes_dir)/ON/ON-NEGF-share


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

all-xk-gnu: rmg-xk-gnu

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

rmg-xk-gnu: RMG/Headers/make_conf.h
	@echo "#define MPI 1" > RMG/Headers/arch.h
	@echo "#define LINUX 1" >> RMG/Headers/arch.h
	cd lib/libxc/; $(MAKE) -f Make.rmg
	@cd RMG; $(clean-global)
	@echo "#define BUILD_DATE \"`date \"+%c\"`\"" >> RMG/Headers/arch.h
	cd RMG; $(MAKE) -f Make.xk-gnu 2>&1 | tee .build.log

rmg-aix: RMG/Headers/make_conf.h
	@echo "#define MPI 1" > RMG/Headers/arch.h
	@echo '#define AIX 1' >> RMG/Headers/arch.h
	cd lib/libxc/; $(MAKE) -j 8 -f Make.rmg
	@cd RMG; $(clean-global)
	@echo "#define BUILD_DATE \"`date \"+%c\"`\"" >> RMG/Headers/arch.h
	cd RMG;gmake -f Make.aix 2>&1 | tee .build.log


# Order-N targets go here
on-linux: 
	@echo "#define LINUX 1" > $(ON_dir)/Headers/arch.h
	@echo "#define MPI 1" >> $(ON_dir)/Headers/arch.h
	@echo "#define HYBRID_MODEL 0" >> $(ON_dir)/Headers/arch.h;
	@cd $(ON_dir); $(clean-global); $(clean-on-negf-share)
	cd $(ON_dir); $(MAKE) -j 8 -f Make.linux 2>&1 | tee .build.log

on-xt: 
	@echo "#define LINUX 1" > $(ON_dir)/Headers/arch.h
	@echo "#define MPI 1" >> $(ON_dir)/Headers/arch.h
	@echo "#define HYBRID_MODEL 0" >> $(ON_dir)/Headers/arch.h;
	@cd $(ON_dir); $(clean-global); $(clean-on-negf-share)
	cd $(ON_dir); $(MAKE) -f Make.xt 2>&1 | tee .build.log


on-aix: 
	@echo '#define AIX_MPI 1' > $(ON_dir)/Headers/arch.h
	@echo "#define PARALLEL_MESSAGE 1" >> $(ON_dir)/Headers/arch.h
	@echo "#define HYBRID_MODEL 0" >> $(ON_dir)/Headers/arch.h;
	@cd $(ON_dir); $(clean-global); $(clean-on-negf-share)
	cd $(ON_dir); gmake -f Make.aix 2>&1 | tee .build.log


# Targets for NEGF compilation go here
NEGF-linux: 
	@echo '#define LINUX 1' > $(NEGF_dir)/Headers/arch.h
	@echo "#define MPI 1" >> $(NEGF_dir)/Headers/arch.h
	@echo "#define HYBRID_MODEL 0" >> $(NEGF_dir)/Headers/arch.h;
	@echo "#define GPU_ENABLE 0" >> $(NEGF_dir)/Headers/arch.h;
	@cd $(NEGF_dir); $(clean-global); $(clean-on-negf-share)
	cd $(NEGF_dir); $(MAKE) -j 8 -f Make.linux 2>&1 | tee .build.log

NEGF-xt: 
	@echo '#define LINUX 1' > $(NEGF_dir)/Headers/arch.h
	@echo "#define MPI 1" >> $(NEGF_dir)/Headers/arch.h
	@echo "#define HYBRID_MODEL 0" >> $(NEGF_dir)/Headers/arch.h;
	@echo "#define GPU_ENABLE 1" >> $(NEGF_dir)/Headers/arch.h;
	@cd $(NEGF_dir); $(clean-global); $(clean-on-negf-share)
	cd $(NEGF_dir); $(MAKE) -f Make.xt 2>&1 | tee .build.log

NEGF-aix: 
	@echo '#define AIX_MPI 1' > $(NEGF_dir)/Headers/arch.h
	@echo "#define PARALLEL_MESSAGE 1" >> $(NEGF_dir)/Headers/arch.h
	@echo "#define HYBRID_MODEL 0" >> $(NEGF_dir)/Headers/arch.h;
	@echo "#define GPU_ENABLE 0" >> $(NEGF_dir)/Headers/arch.h;
	@cd $(NEGF_dir); $(clean-global); $(clean-on-negf-share)
	cd $(NEGF_dir); gmake -f Make.aix 2>&1 | tee .build.log

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

clean-all: clean-common clean-rmg clean-on clean-NEGF clean-libxc

clean-libxc:
	cd lib/libxc; make clean 



RMG/Headers/make_conf.h:
	@echo "/****** Compile-time options for real-space code are set here  ******/" > RMG/Headers/make_conf.h;
	@echo "" >> RMG/Headers/make_conf.h;
	@echo "/* To enable MPI/PThreads hybrid model, experimental for now */ " >> RMG/Headers/make_conf.h;
	@echo "#define HYBRID_MODEL 0" >> RMG/Headers/make_conf.h;
	@echo "/* Number of threads in MPI/PThreads mode */ " >> RMG/Headers/make_conf.h;
	@echo "" >> RMG/Headers/make_conf.h;
	@echo "/* Gamma point only set to 1, otherwise, 0 */" >> RMG/Headers/make_conf.h;
	@echo "#define GAMMA_PT 1" >> RMG/Headers/make_conf.h;
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
	@echo "/* Use Magma libs. */" >> RMG/Headers/make_conf.h;
	@echo "#define MAGMA_LIBS 0" >> RMG/Headers/make_conf.h;
	@echo "";
	@echo "/* Use Async trade_images. */" >> RMG/Headers/make_conf.h;
	@echo "#define ASYNC_TRADES 0" >> RMG/Headers/make_conf.h;
	@echo "";
	@echo "/* Number of OMP threads to use (internal only). If 0 then none will be used. */" >> RMG/Headers/make_conf.h;
	@echo "#define RMG_OMP_THREADS 0" >> RMG/Headers/make_conf.h;
	@echo "";
	@echo "ERROR: File RMG/Headers/make_conf.h does not exist"; 
	@echo "Headers/make_conf.h is set to default, configure it for your system before compiling !!! ";
	@echo "";
	exit 1; 

