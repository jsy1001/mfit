# $Id: makefile,v 1.12 2003/07/18 17:50:23 jsy1001 Exp $
#
# Makefile for building mfit

SHELL = /bin/sh

# Sun Workshop Fortran 95
#F90 = /opt/SUNWspro/bin/f90
#FFLAGC = -g -dalign
#FFLAGL = -dalign -lf77compat
#pgplot_libs = -lpgplot -lX11
#fitsio_libs = -lfitsio
#pda_libs = -L/star/lib -lpda -lemsf -lems -lcnf

# NAGWare Fortran 95 on Solaris
F90 = /opt/NAGWare/bin/f95
FFLAGC = -g -mismatch
FFLAGL = -lF77 -lM77
pgplot_libs = -lpgplot -lX11
fitsio_libs = -lfitsio
pda_libs = -L/star/lib -lpda -lemsf -lems -lcnf

# NAGWare Fortran 95 on Linux
#F90 = /usr/local/bin/f95
#FFLAGC = -g -mismatch
# Not sure about next line
#FFLAGL = -lg2c -lm
#pgplot_libs = -lpgplot -lX11
#fitsio_libs = -lfitsio
#pda_libs = -L/star/lib -lpda -lemsf -lems -lcnf

# *** Modify above this line for your system ***

OBJECTS = maths.o fit.o visibility.o inout.o plot.o model.o \
	gamma.o rjbesl.o
MODULES = maths.mod fit.mod visibility.mod inout.mod plot.mod model.mod


all: mfit clfit calc ;

mfit: main.o $(OBJECTS)
	$(F90) $(FFLAGL) $^ -o $@ $(pgplot_libs) $(fitsio_libs) $(pda_libs)

clfit: clfit.o f2kcli.o $(OBJECTS)
	$(F90) $(FFLAGL) $^ -o $@ $(pgplot_libs) $(fitsio_libs) $(pda_libs)

calc: calc.o maths.o gamma.o rjbesl.o
	$(F90) $(FFLAGL) $^ -o $@ $(pda_libs)

clean:
	rm -f *.o *.mod *.lst

# source files containing module definitions must be compiled before source
# files that USE those modules
main.o: main.f90 inout.mod plot.mod visibility.mod fit.mod model.mod
	$(F90) -c $(FFLAGC) main.f90

clfit.o: clfit.f90 f2kcli.mod inout.mod plot.mod visibility.mod fit.mod model.mod
	$(F90) -c $(FFLAGC) clfit.f90

calc.o: calc.f90 maths.mod
	$(F90) -c $(FFLAGC) calc.f90

f2kcli.o f2kcli.mod: f2kcli.f90
	$(F90) -c $(FFLAGC) f2kcli.f90

maths.o maths.mod: maths.f90
	$(F90) -c $(FFLAGC) maths.f90

fit.o fit.mod: fit.f90 maths.mod visibility.mod model.mod
	$(F90) -c $(FFLAGC) fit.f90

visibility.o visibility.mod: visibility.f90 maths.mod
	$(F90) -c $(FFLAGC) visibility.f90

inout.o inout.mod: inout.f90
	$(F90) -c $(FFLAGC) inout.f90

plot.o plot.mod: plot.f90 model.mod fit.mod
	$(F90) -c $(FFLAGC) plot.f90

model.o model.mod: model.f90
	$(F90) -c $(FFLAGC) model.f90

gamma.o: gamma.f
	$(F90) -c $(FFLAGC) gamma.f

rjbesl.o: rjbesl.f
	$(F90) -c $(FFLAGC) rjbesl.f
