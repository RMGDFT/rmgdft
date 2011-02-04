#************************** SVN Revision Information **************************
#**    $Id: Makefile 1242 2011-02-02 18:55:23Z luw $    **
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


aix: 
	@echo '#define AIX_MPI 1' > Headers/arch.h
	@echo "#define PARALLEL_MESSAGE 1" >> Headers/arch.h
	gmake -f Make.aix

sunfire: Headers/make_conf.h
	@echo '#define LINUX 1' > Headers/arch.h
	@echo "#define MPI 1" >> Headers/arch.h
	$(MAKE) -f Make.sunfire

linux: 
	@echo '#define LINUX 1' > Headers/arch.h
	@echo "#define MPI 1" >> Headers/arch.h
	$(MAKE) -f Make.linux

linux_mpich:
	@echo '#define LINUX 1' > Headers/arch.h
	@echo "#define MPI 1" >> Headers/arch.h
	$(MAKE) -f Make.linux_mpich

xt: 
	@echo '#define LINUX 1' > Headers/arch.h
	@echo "#define MPI 1" >> Headers/arch.h
	$(MAKE) -f Make.xt


clean:
	find . \( -name '*.o' -o -name '*.oo' \) -exec rm {} \;
	find ../Common \( -name '*.o' -o -name '*.oo' \) -exec rm {} \;
	rm -f SvnRev/svnrev Headers/svnrev.h Headers/arch.h Headers/grid.h Make.vars loadmap core mppcore rmg-on

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
