#### May need to edit this section ####

# These may be overridden using "debug=" and "release=" command-line arguments
debug = 1
release = 0

includePath = Split('/opt/local/include /star/include /sw/include')
libPath = Split('/opt/local/lib /star/lib /sw/lib')

# Path to Sun F95 (if installed); used to auto-set compiler switches
sun_f95 = '/opt/SUNWspro/bin/f95'
# Path to NAGWare F95 (if installed); used to auto-set compiler switches
nagw_f95 = '/opt/NAGWare/bin/f95'

# Fortran 90/95 compilers to search for, in preferred order
f95List = ['g95', sun_f95, nagw_f95, 'f95', 'f90']

####### End of editable section #######


env = Environment()

# Parse "debug=" and "release=" command-line arguments (defaults used if not
# specified are in configuration section above)
debug = int(ARGUMENTS.get('debug', debug))
release = int(ARGUMENTS.get('release', release))

# If SunOS, modify PATH so Sun compilers are found
# Not necessary if full f95 path supplied above
if Platform() == 'sunos':
    path = ['/opt/SUNWspro/bin', '/usr/local/bin',
            '/bin', '/usr/bin', '/usr/ccs/bin']
    env.Append(ENV={'PATH' : path})

# Find fortran compiler and set command-line switches to use for it
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
    env.Append(FORTRANFLAGS=['-g'])
if f95 == 'g95':
    env.Append(FORTRANFLAGS=['-fno-second-underscore'])
    f2kcli = 'f2kcli.f90'
elif f95 == sun_f95:
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
    env.Append(FORTRANFLAGS=['-mismatch'])
    # Prevent 'undefined symbol' errors when linking to libraries built with
    # Sun fortran compiler
    env.Append(LIBPATH=['/opt/SUNWspro/lib'])
    env.Append(LIBS=['F77', 'M77'])
    # need different f2kcli for this compiler
    f2kcli = 'f2kcli_nagw.f90'
else:
    f2kcli = 'f2kcli.f90'
    
# Prepend to include/library paths from variables above
env.Append(FORTRANPATH=includePath)
# FORTRANPATH seems not to be used for .f90
env.Prepend(F90PATH=includePath)
env.Prepend(LIBPATH=libPath)
if Platform() != 'win32':
    libUnixPath = ':'.join(libPath)
    env.Append(ENV={'LD_RUN_PATH' : libUnixPath})

# Configure libraries to link in
conf = Configure(env)
pdaLibs = ['pda', 'emsf', 'ems', 'cnf']
if conf.CheckLib('starmem'):
    pdaLibs += ['starmem']
pgplotLibs = ['pgplot', 'X11']
fitsioLibs = ['cfitsio']
if conf.CheckLib('socket'):
    fitsioLibs += ['socket', 'nsl']
fftwLibs = ['fftw3']
if Platform() == 'sunos':
    fftwLibs += ['m']

mfitLibs = env['LIBS'] + pdaLibs + pgplotLibs + fitsioLibs + fftwLibs
if conf.CheckLib('nag'):
    conf.env.Append(CPPDEFINES='HAVE_NAG=1')
    mfitLibs += ['nag']
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
sources['mplot'] = ['modelplot.f90',
                    'maths.f90', 'model.f90',
                    'gamma.f', 'rjbesl.f', 'fitsimage.f90', 'inout.f90']
sources['calc'] = ['calc.f90',
                   'maths.f90', 'gamma.f', 'rjbesl.f']
libs = {}
libs['mfit'] = mfitLibs
libs['clfit'] = mfitLibs
libs['mplot'] = env['LIBS'] + pdaLibs + pgplotLibs + fitsioLibs + fftwLibs
libs['calc'] = env['LIBS'] + pdaLibs
objects = {}
# ...object files
for key in sources.keys():
    allobjs = env.Object(sources[key])
    # filter out mod files
    objects[key] = filter(lambda o: str(o)[-4:] != '.mod', allobjs)

# ...executables
for key in objects.keys():
    env.Program(key, objects[key], LIBS=libs[key])

# Local Variables:
# mode: python
# End:
