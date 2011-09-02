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


# Order-N targets go here
on-linux: 
	@echo "#define LINUX 1" > Headers/arch.h
	@echo "#define MPI 1" >> Headers/arch.h
	@echo "#define HYBRID_MODEL 0" >> Headers/arch.h
	@echo "#define THREADS_PER_NODE 1" >> Headers/arch.h
	@$(clean-global); $(clean-on-negf-share)
	$(MAKE) -j 8 -f Make.linux 2>&1 | tee .build.log

on-xt: 
	@echo "#define LINUX 1" > Headers/arch.h
	@echo "#define MPI 1" >> Headers/arch.h
	@echo "#define HYBRID_MODEL 0" >> Headers/arch.h
	@echo "#define THREADS_PER_NODE 1" >> Headers/arch.h
	@$(clean-global); $(clean-on-negf-share)
	$(MAKE) -f Make.xt 2>&1 | tee .build.log

on-aix: 
	@echo '#define AIX_MPI 1' > Headers/arch.h
	@echo "#define PARALLEL_MESSAGE 1" >> Headers/arch.h
	@echo "#define HYBRID_MODEL 0" >> Headers/arch.h
	@echo "#define THREADS_PER_NODE 1" >> Headers/arch.h
	@$(clean-global); $(clean-on-negf-share)
	gmake -f Make.aix 2>&1 | tee .build.log


#Clean targets

.PHONY: clean-common
clean-common:
	find $(GLOBAL_MODULES)  \( -name '*.o' -o -name '*.oo' \) -exec rm {} \;

clean: clean-common
	find .  \( -name '*.o' -o -name '*.oo' \) -exec rm {} \;

clean-libxc:
	cd ../lib/libxc; make clean 

