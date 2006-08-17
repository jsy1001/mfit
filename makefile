# $Id: makefile,v 1.20 2006/08/17 13:34:32 jsy1001 Exp $
#
# Makefile for building mfit

SHELL = /bin/sh

ifeq ($(OSTYPE),solaris)
  # Sun Workshop Fortran 95
  F90 = /opt/SUNWspro/bin/f90
  FLINK = /opt/SUNWspro/bin/f90
  FFLAGC = -C -g -dalign -I/opt/local/include
  FFLAGL = -dalign -lf77compat
  #FFLAGC = -g -fast -dalign -I/opt/local/include
  #FFLAGL = -fast -dalign -lf77compat
  f2kcli_src = f2kcli.f90

  # NAGWare Fortran 95 on Solaris
  #F90 = /opt/NAGWare/bin/f95
  #FLINK = /opt/NAGWare/bin/f95
  #FFLAGC = -g -mismatch -I/opt/local/include
  #FFLAGL = -lF77 -lM77
  #f2kcli_src = f2kcli_nagw.f90

  # Libraries needed on Solaris
  pgplot_libs = -lpgplot -lX11
  fitsio_libs = -lfitsio
  pda_libs = -L/star/lib -lpda -lemsf -lems -lcnf
  nag_libs = -lnag
  fftw_libs = -lfftw3 -lm
else
  # NAGWare Fortran 95 on Linux
  F90 = /usr/local/bin/f95
  FLINK = /usr/local/bin/f95
  FFLAGC = -g -mismatch -I/usr/local/include
  #FFLAGL = -lg2c -lm
  # uncomment next line to build semi-static binary (for systems without f95):
  FFLAGL = -unsharedf95 -lg2c -lm
  f2kcli_src = f2kcli_nagw.f90

  # Libraries needed on Linux
  pgplot_libs = -L/usr/X11R6/lib -L/usr/local/lib -lpgplot -lX11
  fitsio_libs = -L/star/lib -lfitsio
  pda_libs = -L/star/lib -lpda -lemsf -lems -lcnf
  nag_libs = -lnag
  fftw_libs = -L/usr/local/lib -lfftw3 -lm
endif

# *** Modify above this line for your system ***

OBJECTS = maths.o fit.o visibility.o inout.o plot.o model.o \
  gamma.o rjbesl.o marginalise.o
MODULES = maths.mod fit.mod visibility.mod inout.mod plot.mod model.mod marginalise.mod


all: mfit clfit calc mplot ;

mfit: main.o $(OBJECTS)
	$(FLINK) $(FFLAGL) $^ -o $@ $(pgplot_libs) $(fitsio_libs)  $(nag_libs) $(pda_libs) $(fftw_libs)

clfit: clfit.o f2kcli.o $(OBJECTS)
	$(FLINK) $(FFLAGL) $^ -o $@ $(pgplot_libs) $(fitsio_libs)  $(nag_libs) $(pda_libs) $(fftw_libs)

calc: calc.o maths.o gamma.o rjbesl.o
	$(FLINK) $(FFLAGL) $^ -o $@ $(pda_libs)

mplot: modelplot.f90 model.o maths.o fitsimage.o gamma.o rjbesl.o inout.o
	$(FLINK) $(FFLAGL) $^ -o $@ $(pgplot_libs) $(fitsio_libs) $(pda_libs) $(fftw_libs)

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

f2kcli.o f2kcli.mod: $(f2kcli_src)
	$(F90) -c $(FFLAGC) $< -o f2kcli.o

maths.o maths.mod: maths.f90
	$(F90) -c $(FFLAGC) maths.f90

fit.o fit.mod: fit.f90 maths.mod visibility.mod model.mod
	$(F90) -c $(FFLAGC) fit.f90

visibility.o visibility.mod: visibility.f90 maths.mod model.mod
	$(F90) -c $(FFLAGC) visibility.f90

inout.o inout.mod: inout.f90
	$(F90) -c $(FFLAGC) inout.f90

plot.o plot.mod: plot.f90 model.mod fit.mod marginalise.mod visibility.mod
	$(F90) -c $(FFLAGC) plot.f90

model.o model.mod: model.f90 maths.mod
	$(F90) -c $(FFLAGC) model.f90

marginalise.o marginalise.mod: marginalise.f90 fit.mod model.mod
	$(F90) -c $(FFLAGC) marginalise.f90

gamma.o: gamma.f
	$(F90) -c $(FFLAGC) gamma.f

rjbesl.o: rjbesl.f
	$(F90) -c $(FFLAGC) rjbesl.f

modelplot.o: modelplot.f90 model.mod maths.mod fitsimage.mod inout.mod 
	$(F90) -c $(FFLAGC) modelplot.f90

fitsimage.o: fitsimage.f90
	$(F90) -c $(FFLAGC) fitsimage.f90
