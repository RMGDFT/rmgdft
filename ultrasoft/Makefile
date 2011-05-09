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

default:
	@echo "Type one of the following:"
	@echo "  make aix                  for 64 bit AIX using MPI"
	@echo "  make sunfire              for sunfire opteron using MPI"
	@echo "  make linux                for linux using MPI"
	@echo "  make xt                   for CRAY XT using MPI"
	@echo " "


aix: Headers/make_conf.h
	@echo "#define MPI 1" > Headers/arch.h
	@echo '#define AIX 1' >> Headers/arch.h
	gmake -f Make.aix

sunfire: Headers/make_conf.h
	@echo "#define MPI 1" > Headers/arch.h
	@echo '#define LINUX 1' >> Headers/arch.h
	cd ../lib/libxc/; $(MAKE) -f Make.rmg
	$(MAKE) -f Make.sunfire

linux: Headers/make_conf.h
	@echo "#define MPI 1" > Headers/arch.h
	@echo '#define LINUX 1' >> Headers/arch.h
	cd ../lib/libxc/; $(MAKE) -j 8 -f Make.rmg
	$(MAKE) -j 8 -f Make.linux

linux_mpich: Headers/make_conf.h
	@echo "#define MPI 1" > Headers/arch.h
	@echo '#define LINUX 1' >> Headers/arch.h
	cd ../lib/libxc/; $(MAKE) -f Make.rmg
	$(MAKE) -j 8 -f Make.linux_mpich

xt: Headers/make_conf.h
	@echo "#define MPI 1" > Headers/arch.h
	@echo '#define LINUX 1' >> Headers/arch.h
	cd ../lib/libxc/; $(MAKE) -f Make.rmg
	$(MAKE) -f Make.xt

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
	@echo "/* How many times the fine grid is finer than the coarse grid " >>Headers/make_conf.h;
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


clean:
	find . \( -name '*.o' -o -name '*.oo' \) -exec rm {} \;
	find ../Common \( -name '*.o' -o -name '*.oo' \) -exec rm {} \;
	rm -f SvnRev/svnrev Headers/svnrev.h Headers/arch.h Headers/grid.h Make.vars loadmap core mppcore rmg

tar:
	rm -f mg.tar.*
	tar -cf mg.tar CHANGELOG README* LICENSE *.h *.c *.f Make.* *.dxx doc/*
	gzip mg.tar

alltar:
	rm -f mg.tar.*
	tar -cf mg.tar CHANGELOG README* LICENSE *.h *.c *.f Makefile Make.* pseudo.* in.* *.dxx doc/*
	gzip mg.tar                                                                                               

htmldoc:
	@echo "#define M_LINUX_MPI 1" > arch.h
	@echo "#define NPES "$(NPES) >> arch.h
	doc++ -v -eo 12pt --header docheader.dxx --footer docfooter.dxx -a -S -m -j --dir doc doc.dxx
DOC:
	doc.bat
