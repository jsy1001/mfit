project('mfit', 'fortran')

cfitsio_dep = dependency('cfitsio')
fftw_dep = dependency('fftw3')

mfit_lib = static_library('mfit', 'maths.f90', 'search.f90',
                          'visibility.f90', 'inout.f90', 'model.f90',
                          'gamma.f', 'rjbesl.f', 'maths_pda.f',
                          'wrap.f90', 'bayes.f90',
                          fortran_args : ['-I/usr/include'],
                          dependencies : [cfitsio_dep, fftw_dep])

executable('mfit', 'main.f90',
           'fit.f90', 'plot.f90', 'postplot.f90', 'marginalise.F90',
           'fit_pda.f', 'pda_xermsg.f', 'gmst.f', 'dranrm.f',
           link_with : mfit_lib,
           link_args : ['-lpgplot'],
           install : true)

executable('clfit', 'clfit.f90',
           'fit.f90', 'plot.f90', 'postplot.f90', 'marginalise.F90',
           'fit_pda.f', 'pda_xermsg.f', 'gmst.f', 'dranrm.f',
           'f2kcli.f90',
           link_with : mfit_lib,
           link_args : ['-lpgplot'],
           install : true)

executable('mplot', 'modelplot.f90',
           'fitsimage.f90',
           link_with : mfit_lib,
           link_args : ['-lpgplot'],
           install : true)

executable('calc', 'calc.f90',
           'maths.f90', 'gamma.f', 'rjbesl.f', 'maths_pda.f')

executable('clnest', 'clnest.f90',
           'nestwrap.f90', 'readmc.f90', 'f2kcli.f90',
           link_with : mfit_lib,
           build_by_default : false)

executable('binnest', 'binnest.f90',
           'nestwrap.f90', 'readmc.f90', 'f2kcli.f90',
           link_with : mfit_lib,
           build_by_default : false)

