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



# RMG, real-space US global grid code:
rmg-linux: Headers/make_conf.h
	@echo "#define MPI 1" > Headers/arch.h
	@echo "#define LINUX 1" >> Headers/arch.h
	cd ../lib/libxc/; $(MAKE) -f Make.rmg
	@$(clean-global)
	@echo "#define BUILD_DATE \"`date \"+%c\"`\"" >> Headers/arch.h
	$(MAKE) -j 8 -f Make.linux 2>&1 | tee .build.log
	

rmg-xt: Headers/make_conf.h
	@echo "#define MPI 1" > Headers/arch.h
	@echo "#define LINUX 1" >> Headers/arch.h
	cd ../lib/libxc/; $(MAKE) -f Make.rmg
	@$(clean-global)
	@echo "#define BUILD_DATE \"`date \"+%c\"`\"" >> Headers/arch.h
	$(MAKE) -f Make.xt 2>&1 | tee .build.log

rmg-aix: Headers/make_conf.h
	@echo "#define MPI 1" > Headers/arch.h
	@echo "#define LINUX 1" >> Headers/arch.h
	cd ../lib/libxc/; $(MAKE) -f Make.rmg
	@$(clean-global)
	@echo "#define BUILD_DATE \"`date \"+%c\"`\"" >> Headers/arch.h
	gmake -f Make.aix 2>&1 | tee .build.log


#Clean targets

.PHONY: clean-common
clean-common:
	find $(GLOBAL_MODULES)  \( -name '*.o' -o -name '*.oo' \) -exec rm {} \;

clean: clean-common
	find .  \( -name '*.o' -o -name '*.oo' \) -exec rm {} \;

clean-libxc:
	cd ..lib/libxc; make clean 



Headers/make_conf.h:
	@echo "/****** Compile-time options for real-space code are set here  ******/" > Headers/make_conf.h;
	@echo "" >> Headers/make_conf.h;
	@echo "/* Number of processors */" >> Headers/make_conf.h;
	@echo "#define NPES 2" >> Headers/make_conf.h;
	@echo "" >> Headers/make_conf.h;
	@echo "/* 3D processor grid PE_X*PE_Y*PE_Z must equal NPES */ " >> Headers/make_conf.h;
	@echo "#define PE_X 1" >> Headers/make_conf.h;
	@echo "#define PE_Y 1" >> Headers/make_conf.h;
	@echo "#define PE_Z 2" >> Headers/make_conf.h;
	@echo "" >> Headers/make_conf.h;
	@echo "/* Gamma point only set to 1, otherwise, 0 */" >> Headers/make_conf.h;
	@echo "#define GAMMA_PT 1" >> Headers/make_conf.h;
	@echo "" >> Headers/make_conf.h;
	@echo "/* Number of points in coarse grid in x,y and z directions*/" >> Headers/make_conf.h;
	@echo "#define NX_GRID 32" >> Headers/make_conf.h;
	@echo "#define NY_GRID 32" >> Headers/make_conf.h;
	@echo "#define NZ_GRID 32" >> Headers/make_conf.h;
	@echo "" >> Headers/make_conf.h;
	@echo "/* How many times the fine grid is finer than the coarse grid " >> Headers/make_conf.h;
	@echo " * All three numbers have to be the same */" >> Headers/make_conf.h;
	@echo "#define FG_NX 2" >> Headers/make_conf.h;
	@echo "#define FG_NY 2" >> Headers/make_conf.h;
	@echo "#define FG_NZ 2" >> Headers/make_conf.h;
	@echo "" >> Headers/make_conf.h;
	@echo "/* Set this to 0 to turn off memory Smart-ALLOCation. (salloc.c, salloc.h)" >> Headers/make_conf.h;
	@echo " * (there is no significant time improvement)*/" >> Headers/make_conf.h;
	@echo "#define USE_SALLOC 1" >> Headers/make_conf.h;
	@echo "" >> Headers/make_conf.h;
	@echo "/* Set this to 1 if you want to use finite difference method for calculating" >> Headers/make_conf.h;
	@echo " * derivatives of beta. This is faster since it avoids doing 3 backwards fourier" >> Headers/make_conf.h;
	@echo " * transforms per ion, but it may not be very accurate since the finite diff" >> Headers/make_conf.h;
	@echo " * derivative is done on the coarse grid." >> Headers/make_conf.h;
	@echo " * Leave this set to 0 unless you know what you are doing */" >> Headers/make_conf.h;
	@echo "#define FDIFF_BETA 0" >> Headers/make_conf.h;
	@echo "" >> Headers/make_conf.h;
	@echo "/* Extra fine timers, may cause some slowdown, but they are generally useful */" >> Headers/make_conf.h;
	@echo "#define MD_TIMERS 1" >> Headers/make_conf.h;
	@echo "";
	@echo "ERROR: File Headers/make_conf.h does not exist"; 
	@echo "Headers/make_conf.h is set to default, configure it for your system before compiling !!! ";
	@echo "";
	exit 1; 

