import sys
from numpy.distutils.core import Extension

# :TODO: determine correct flags for include 'fftw3.f'
sys.argv.extend(["config_fc", "--f90flags='-I/usr/include'"])

ext1 = Extension(name='mfit',
                 sources=['model.f90', 'visibility.f90', 'bayes.f90',
                          'fit.f90', 'inout.f90'],
                 libraries=['mfit', 'cfitsio', 'fftw3'])

# note module order matters
nowrap_sources = ['maths.f90', 'search.f90', 'component.f90', 'wrap.f90',
                  'gamma.f', 'rjbesl.f', 'maths_pda.f', 'fit_pda.f',
                  'pda_xermsg.f', 'gmst.f', 'dranrm.f']

if __name__ == "__main__":
    from numpy.distutils.core import setup
    setup(name = 'mfit_py',
          description       = "MFIT Python Interface",
          author            = "John Young",
          author_email      = "jsy1001@cam.ac.uk",
          libraries = [('mfit', dict(sources=nowrap_sources))],
          ext_modules = [ext1]
          )
