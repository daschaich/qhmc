require 'common'
require 'gaugeact'

--trace(true)
latsize = { 4, 4, 4, 8 }
mass = 0.0008
prec = 1
restart = 500
resid = 1e-12
opts = { prec=prec, restart=restart }

--qopqdp.defaultNc(4); fn = nil
Nc = qopqdp.defaultNc()
Ns = 4

L = qopqdp.lattice(latsize)
qopqdp.profile(profile or 0)
qopqdp.verbosity(0)

seed = seed or os.time()
qopqdp.seed(seed)

function getplaq(g)
  local ps,pt = g:action({plaq=1})
  local nd = #qopqdp.lattice()
  local norm_t = (nd - 1) * g:lattice():volume()
  local norm_s = 0.5 * (nd - 2) * norm_t
  ps = ps / norm_s
  pt = pt / norm_t
  printf("plaq ss: %-8g  st: %-8g  tot: %-8g\n", ps, pt, 0.5*(ps+pt))
end

g = qopqdp.gauge()
if fn then g:load(fn)
else
  --g:unit()
  g:random()
end

getplaq(g)
-- coulomb(j_decay, error, max iterations, overrelaxation param)
-- note that here 0->x, ..., 3->t
--g:coulomb(3, 1e-7, 1000, 1.2)
--getplaq(g)

h = qopqdp.hisq()
--h:printcoeffs()
h:set(g, 1, 0, prec)

local pt = {0,0,0,0}
printf("src point: %i %i %i %i\n", pt[1], pt[2], pt[3], pt[4])
src = h:quark()
dest = h:quark()
src:point(pt,1,1)
dest:zero()

t0 = qopqdp.dtime()
h:solve({dest}, src, {mass}, resid, "all", opts)
dt = qopqdp.dtime() - t0
mf = 1e-6 * h:flops() / dt
printf("its: %g  secs: %g  Mflops: %g\n", h:its(), dt, mf)
