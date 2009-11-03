#### May need to edit this section ####

# These may be overridden using "debug=" and "release=" command-line arguments
debug = 1
release = 0

# Locations of 64-bit files appear first,
# since they don't exist on 32-bit machines
includePath = Split('/opt/local/include /soft/star64/include /opt/star/include /star/include /sw/include /usr/include ./MultiNest_v2.7')
libPath = Split('/opt/local/lib /usr/lib64 /soft/star64/lib /opt/star/lib /star/lib /sw/lib /usr/lib /usr/X11R6/lib64 /usr/X11R6/lib ./MultiNest_v2.7')

# Path to Sun F95 (if installed); used to auto-set compiler switches
sun_f95 = '/opt/SUNWspro/bin/f95'
# Path to NAGWare F95 (if installed); used to auto-set compiler switches
nagw_f95 = '/usr/local/nag/bin_amd64/f95'
#nagw_f95 = '/opt/NAGWare/bin/f95'

# Fortran 90/95 compilers to search for, in preferred order
f95List = ['g95', sun_f95, nagw_f95, 'f95', 'f90']

####### End of editable section #######


# Inherit complete environment, including PATH
import os
env = Environment(ENV = os.environ)
conf = Configure(env)

# Parse command-line arguments
debug = int(ARGUMENTS.get('debug', debug))
release = int(ARGUMENTS.get('release', release))
user_f95 = ARGUMENTS.get('f95', None)

# Find fortran compiler and set command-line switches to use for it
if user_f95 is not None:
    f95 = user_f95
else:    
    f95 = env.Detect(f95List)
    if f95 is None:
        print "No f90/f95 compiler found in %s" % env['ENV']['PATH']
        Exit(1)
print "Building on '%s' using %s" % (Platform(), f95)
print "debug=%s  release=%s" % (debug, release)
env.Replace(FORTRAN=f95, LINK=f95)
env.Replace(FORTRANFLAGS=[])
env.Replace(LINKFLAGS=[])
if debug:
    # most compilers understand -g
    env.Append(FORTRANFLAGS=['-g'])
if f95 == 'g95':
    print "Configuring for G95"
    env.Append(FORTRANFLAGS=['-fno-second-underscore'])
    f2kcli = 'f2kcli.f90'
elif f95 == sun_f95:
    print "Configuring for Sun WorkShop Fortran 95"
    env.Append(FORTRANFLAGS=['-dalign'])
    if release:
        env.Append(FORTRANFLAGS=['-fast'])
        env.Append(LINKFLAGS=['-fast'])
    if debug:
        env.Append(FORTRANFLAGS=['-C'])
    env.Append(LINKFLAGS=['-dalign'])
    env.Append(LIBS=['f77compat'])
    f2kcli = 'f2kcli.f90'
elif f95 == nagw_f95:
    print "Configuring for NAGWare Fortran 95"
    env.Append(FORTRANFLAGS=['-mismatch'])
    if release:
        env.Append(FORTRANFLAGS=['-O'])
        env.Append(LINKFLAGS=['-O'])
    if debug:
        env.Append(FORTRANFLAGS=['-C'])
    # Prevent 'undefined symbol' errors when linking to libraries built with
    # Sun fortran compiler
    env.Append(LIBPATH=['/opt/SUNWspro/lib'])
    if conf.CheckLib('F77'):
        # NB CheckLib() appends to env['LIBS'] if found
        env.Append(LIBS=['M77']) # need both F77 and M77
    # need different f2kcli for this compiler
    f2kcli = 'f2kcli_nagw.f90'
    # Ignores LD_RUN_PATH, so pass "-rpath <path>" to ld
    if env['PLATFORM'] in ['linux','posix']:
        # f95 runs gcc which runs linux ld
        env.Append(LINKFLAGS=['-Wl,-Xlinker', '-Wl,-rpath',
                              '-Wl,-Xlinker', '-Wl,%s' % ':'.join(libPath)])
    elif env['PLATFORM'] == 'sunos':
        # f95 runs cc or gcc (specified at purchase) which runs sunos ld
        # Fortunately cc understands -R directly
        env.Append(LINKFLAGS=['-Wl,-R%s' % ':'.join(libPath)])
else:
    print "Configuring for generic Fortran 9x compiler"
    f2kcli = 'f2kcli.f90'
# Work around inconsistent behaviour between scons versions
env.Replace(F90FLAGS=env['FORTRANFLAGS'])
env.Replace(F90=f95)

# Prepend to include/library paths from variables above
env.Prepend(FORTRANPATH=includePath)
# FORTRANPATH seems not to be used for .f90
env.Prepend(F90PATH=includePath)
env.Prepend(LIBPATH=libPath)
if env['PLATFORM'] != 'win32':
    libUnixPath = ':'.join(env['LIBPATH'])
    env.Append(ENV={'LD_RUN_PATH' : libUnixPath})

# Configure libraries to link in
baseLibs = env.get('LIBS', [])
# include slalib here since its distributed with Starlink
starLibs = ['pda', 'emsf', 'ems', 'cnf', 'sla']
if conf.CheckLib('starmem'):
    # NB CheckLib() appends to env['LIBS'] if found
    starLibs += ['starmem']
pgLibs = ['pgplot', 'X11']
fitsioLibs = ['cfitsio']
if conf.CheckLib('socket'):
    fitsioLibs += ['nsl', 'socket']
fftwLibs = ['fftw3']
if Platform() == 'sunos':
    fftwLibs += ['m']
if conf.CheckLib('nag'):
    conf.env.Append(CPPDEFINES='HAVE_NAG=1')
    nagLibs =  ['nag']
else:
    nagLibs = []
nestLibs = ['nest3', 'lapack']
    
env = conf.Finish()

# Define targets and dependencies...
sources = {}
sources['mfit'] = ['main.f90',
                   'maths.f90', 'fit.f90', 'visibility.f90', 'inout.f90',
                   'plot.f90', 'postplot.f90', 'model.f90',
                   'gamma.f', 'rjbesl.f',
                   'marginalise.F90', 'wrap.f90', 'bayes.f90']
sources['clfit'] = ['clfit.f90',
                   'maths.f90', 'fit.f90', 'visibility.f90', 'inout.f90',
                   'plot.f90', 'postplot.f90', 'model.f90',
                   'gamma.f', 'rjbesl.f',
                   'marginalise.F90', 'wrap.f90', 'bayes.f90'] + [f2kcli]
sources['clnest'] = ['clnest.f90', 'nestwrap.f90', 'readmc.f90',
                   'maths.f90', 'fit.f90', 'visibility.f90', 'inout.f90',
                   'model.f90', 'gamma.f', 'rjbesl.f',
                   'wrap.f90', 'bayes.f90'] + [f2kcli]
sources['binnest'] = ['binnest.f90', 'nestwrap.f90', 'readmc.f90',
                      'wrap.f90', 'model.f90', 'bayes.f90', 'visibility.f90',
                      'maths.f90', 'rjbesl.f'] + [f2kcli]
sources['mplot'] = ['modelplot.f90',
                    'maths.f90', 'model.f90',
                    'gamma.f', 'rjbesl.f', 'fitsimage.f90', 'inout.f90']
sources['calc'] = ['calc.f90',
                   'maths.f90', 'gamma.f', 'rjbesl.f']
libs = {}
libs['mfit'] = baseLibs + starLibs + pgLibs + fitsioLibs + fftwLibs + nagLibs
libs['clfit'] = baseLibs + starLibs + pgLibs + fitsioLibs + fftwLibs + nagLibs
libs['clnest'] = baseLibs + starLibs + fitsioLibs + fftwLibs + nestLibs
libs['binnest'] = baseLibs + nestLibs + starLibs + fftwLibs
libs['mplot'] = baseLibs + starLibs + pgLibs + fitsioLibs + fftwLibs
libs['calc'] = baseLibs + starLibs
objects = {}
# ...object files
for key in sources.keys():
    allobjs = env.Object(sources[key])
    # filter out mod files
    objects[key] = filter(lambda o: str(o)[-4:] != '.mod', allobjs)

# ...executables
for key in objects.keys():
    prog = env.Program(key, objects[key], LIBS=libs[key])
    if key not in ['clnest', 'binnest']:
        Default(prog)

# ...targets for distribution of mfit
import glob
env.Replace(TARFLAGS = '-c -z')
distFiles = Split('README NEWS TODO documentation f2kcli.txt readme.specfun')
distFiles += Split('SConstruct fitgui fitgui_dev')
distFiles += Split('test.oifits test.model')
# don't distribute makefile.coast
for pattern in ['*.f90', '*.F90', '*.f', '*.py']:
    distFiles += glob.glob(pattern)
env.Tar('mfit.tar.gz', distFiles)

# Local Variables:
# mode: python
# End:
