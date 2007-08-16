# $Id: makefile,v 1.22 2007/08/16 16:40:32 jsy1001 Exp $
#
# Makefile for building mfit

SHELL = /bin/sh

# Uncomment the following line if you have the NAg library
HAVE_NAG = true

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
endif
ifeq ($(OSTYPE),darwin8.0)
  # G95 on Mac OS X 10.4
  F90 = g95
  FLINK = g95
  FFLAGC = -g -fno-second-underscore
  FFLAGL = -lM77
  f2kcli_src = f2kcli.f90

  # Libraries needed on Mac OS X
  pgplot_libs = -lpgplot -lX11
  fitsio_libs = -lfitsio
  pda_libs = -L/star/lib -lpda -lemsf -lems -lcnf
  nag_libs = -lnag
  fftw_libs = -lfftw3 -lm
endif
ifeq ($(OSTYPE),linux)
  # NAGWare Fortran 95 on Linux
  F90 = /usr/local/bin/f95
  FLINK = /usr/local/bin/f95
  FFLAGC = -C -g -mismatch -I/usr/local/include
  #FFLAGL = -lg2c -lm
  # uncomment next line to build semi-static binary (for systems without f95):
  FFLAGL = -unsharedf95 -lg2c -lm
  f2kcli_src = f2kcli_nagw.f90

  # Libraries needed on Linux
  pgplot_libs = -L/usr/X11R6/lib -L/usr/local/lib -lpgplot -lX11
  fitsio_libs = -L/usr/lib -lcfitsio
  pda_libs = -L/star/lib -lpda -lemsf -lems -lcnf -lstarmem
  nag_libs = -lnag
  fftw_libs = -L/usr/local/lib -lfftw3 -lm
endif

# *** Modify above this line for your system ***

OBJECTS = maths.o fit.o visibility.o inout.o plot.o postplot.o model.o \
  gamma.o rjbesl.o marginalise.o wrap.o bayes.o
MODULES = maths.mod fit.mod visibility.mod inout.mod plot.mod postplot.mod \
  model.mod marginalise.mod wrap.mod bayes.mod

# Rule to make .o file from .f90 file
%.o %.mod : %.f90
	$(F90) -c $(FFLAGC) $<

# Rule to make .o file from .F90 file (uses preprocessor)
%.o %.mod : %.F90
	$(F90) -fpp -c $(FFLAGC) $(DEFS) $<


all: mfit clfit calc mplot ;

ifeq ($(HAVE_NAG),true)
 mfit_libs = $(pgplot_libs) $(fitsio_libs) $(pda_libs) $(fftw_libs) $(nag_libs)
 DEFS += -DHAVE_NAG
else
 mfit_libs = $(pgplot_libs) $(fitsio_libs) $(pda_libs) $(fftw_libs)
endif

mfit: main.o $(OBJECTS)
	$(FLINK) $(FFLAGL) $^ -o $@ $(mfit_libs)

clfit: clfit.o f2kcli.o $(OBJECTS)
	$(FLINK) $(FFLAGL) $^ -o $@ $(mfit_libs)

calc: calc.o maths.o gamma.o rjbesl.o
	$(FLINK) $(FFLAGL) $^ -o $@ $(pda_libs)

mplot: modelplot.f90 model.o maths.o fitsimage.o gamma.o rjbesl.o inout.o
	$(FLINK) $(FFLAGL) $^ -o $@ \
 $(pgplot_libs) $(fitsio_libs) $(pda_libs) $(fftw_libs)

clean:
	rm -f *.o *.mod *.lst

# source files containing module definitions must be compiled before source
# files that USE those modules
main.o: main.f90 inout.mod model.mod fit.mod plot.mod postplot.mod wrap.mod

clfit.o: clfit.f90 f2kcli.mod inout.mod model.mod fit.mod \
  wrap.mod plot.mod postplot.mod

calc.o: calc.f90 maths.mod

f2kcli.o f2kcli.mod: $(f2kcli_src)
	$(F90) -c $(FFLAGC) $< -o f2kcli.o

wrap.o wrap.mod: wrap.f90

bayes.o bayes.mod: bayes.f90 maths.mod visibility.mod

maths.o maths.mod: maths.f90

fit.o fit.mod: fit.f90 maths.mod bayes.mod wrap.mod model.mod

visibility.o visibility.mod: visibility.f90 maths.mod model.mod

inout.o inout.mod: inout.f90

plot.o plot.mod: plot.f90 maths.mod visibility.mod bayes.mod

postplot.o postplot.mod: postplot.f90 plot.mod bayes.mod wrap.mod \
  model.mod marginalise.mod

model.o model.mod: model.f90 maths.mod

marginalise.o marginalise.mod: marginalise.F90 maths.mod bayes.mod \
  wrap.mod model.mod

modelplot.o: modelplot.f90 model.mod maths.mod fitsimage.mod inout.mod 
