env = Environment()

# :TODO: configuration flags

# G95-specific configuration
env.Append(FORTRAN='g95', LINK='g95')
env.Append(FORTRANFLAGS=['-g', '-fno-second-underscore'])

# Set paths
env.Append(FORTRANPATH=['/star/include', '/sw/include'])
# FORTRANPATH seems not to be used for .f90
env.Append(F90PATH=['/star/include', '/sw/include'])
env.Append(LIBPATH=['/star/lib', '/sw/lib'])
env.Append(ENV={'LD_RUN_PATH' : '/star/lib:/sw/lib'})

##### NO NEED TO EDIT BELOW THIS LINE ####

conf = Configure(env)
pdaLibs = ['pda', 'emsf', 'ems', 'cnf']
if conf.CheckLib('starmem'):
    pdaLibs += ['starmem']
pgplotLibs = ['pgplot', 'X11']

mfitLibs = pdaLibs + pgplotLibs + ['cfitsio', 'fftw3']
if conf.CheckLib('nag'):
    conf.env.Append(CPPDEFINES='HAVE_NAG=1')
    mfitLibs += ['nag']
env = conf.Finish()
sources = ['maths.f90', 'fit.f90', 'visibility.f90', 'inout.f90', 'plot.f90',
           'postplot.f90', 'model.f90', 'gamma.f', 'rjbesl.f',
           'marginalise.F90', 'wrap.f90', 'bayes.f90']
objects = []
for s in sources:
    o = env.Object(s)
    objects += [o[0]] # only want .o, not .mod
env.Program('mfit', ['main.f90']+objects, LIBS=mfitLibs)
env.Program('clfit', ['clfit.f90', 'f2kcli.f90']+objects, LIBS=mfitLibs)
#env.Program('mplot', ['modelplot.f90',
#                      'maths.f90', 'model.f90',
#                      'gamma.f', 'rjbesl.f', 'fitsimage.f90', 'inout.f90'],
#            LIBS=pdaLibs+pgplotLibs+['cfitsio', 'fftw3'])
#env.Program('calc', ['calc.f90', 'maths.f90', 'gamma.f', 'rjbesl.f'],
#            LIBS=pdaLibs)

