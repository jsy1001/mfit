# $Id: makefile,v 1.6 2002/11/18 10:28:58 jsy1001 Exp $
#
# Makefile for mfit on sparc

SHELL = /bin/sh

F90 = /opt/SUNWspro/bin/f90

FFLAGC = -g -C -dalign

package = mfit
prefix = /coast/depot-sparc/$(package)

OBJECTS = main.o maths.o fit.o visibility.o inout.o plot.o model.o \
	gamma.o rjbesl.o
MODULES = maths.mod fit.mod visibility.mod inout.mod plot.mod model.mod

EXES = mfit
TEST_EXES = calc
PACKAGE_DOCS = documentation
REMOVE_TARGETS += $(MODULES)

pda_libs = -lpda -lemsf -lems -lcnf


mfit: $(OBJECTS)
	$(F90) $^ -o $@ `pgplotlink` $(pda_libs) -dalign -lf77compat

calc: calc.o maths.o gamma.o rjbesl.o
	$(F90) $^ -o $@ `$(pda_link)` -dalign -lf77compat

# source files containing module definitions must be compiled before source
# files that USE those modules
main.o: main.f90 inout.mod plot.mod visibility.mod fit.mod model.mod
	$(F90) -c $(FFLAGC) main.f90

maths.mod: maths.f90
	$(F90) -c $(FFLAGC) maths.f90

fit.mod: fit.f90 maths.mod visibility.mod model.mod
	$(F90) -c $(FFLAGC) fit.f90

visibility.mod: visibility.f90 maths.mod
	$(F90) -c $(FFLAGC) visibility.f90

inout.mod: inout.f90
	$(F90) -c $(FFLAGC) inout.f90

plot.mod: plot.f90
	$(F90) -c $(FFLAGC) plot.f90

model.mod: model.f90
	$(F90) -c $(FFLAGC) model.f90

calc.o: calc.f90 maths.mod
	$(F90) -c $(FFLAGC) calc.f90

gamma.o: gamma.f
	$(F90) -c $(FFLAGC) gamma.f

rjbesl.o: rjbesl.f
	$(F90) -c $(FFLAGC) rjbesl.f


include ../build/defaults.mk
