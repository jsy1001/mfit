# $Id: makefile,v 1.1 2002/10/02 16:59:43 jsy1001 Exp $
#
# Makefile for mfit on sparc

SHELL = /bin/sh

F90 = /opt/SUNWspro/bin/f90

FFLAGC = -Xlist -C -g -dalign

package = mfit
prefix = /coast/depot-sparc/$(package)

OBJECTS = main.o maths.o fit.o visibility.o inout.o plot.o
MODULES = maths.mod fit.mod visibility.mod inout.mod plot.mod

EXES = mfit
TEST_EXES = calc
REMOVE_TARGETS += $(MODULES)

mfit: $(OBJECTS)
	$(F90) $^ -o $@ `pgplotlink` -dalign -lnag -lf77compat

calc: calc.o maths.o
	$(F90) $^ -o $@ -dalign -lnag -lf77compat

# source files containing module definitions must be compiled before source
# files that USE those modules
main.o: main.f90 inout.mod plot.mod visibility.mod fit.mod
	$(F90) -c $(FFLAGC) main.f90

maths.mod: maths.f90
	$(F90) -c $(FFLAGC) maths.f90

fit.mod: fit.f90 maths.mod visibility.mod
	$(F90) -c $(FFLAGC) fit.f90

visibility.mod: visibility.f90 maths.mod
	$(F90) -c $(FFLAGC) visibility.f90

inout.mod: inout.f90
	$(F90) -c $(FFLAGC) inout.f90

plot.mod: plot.f90
	$(F90) -c $(FFLAGC) plot.f90

calc.o: calc.f90 maths.mod
	$(F90) -c $(FFLAGC) calc.f90


include ../build/defaults.mk
