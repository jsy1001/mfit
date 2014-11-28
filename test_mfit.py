#!/usr/bin/env/python

import mfit

mfit.inout.read_oi_fits('test.oifits', -1, 0.)
mfit.model.read_model_wrap('test.model')
nvar = mfit.model.model_nvar()
sol, chisqrd, nlposterior, success, info = mfit.fit.minimiser_wrap(nvar)
assert success
print info.rstrip()
