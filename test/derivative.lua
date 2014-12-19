--------------------------------------------------------------------
-- This example is part of the regression test suite.  Functions starting
-- with TEST are part of the test framework and can be ignored.

-- To test chain rule for smeared fields, we start with a gauge field X
-- and some chain field, calculate
--   Re Tr[ f(X) C ]
-- then we add some small value to one element of X, calculate
-- the finite difference and compare to the smearChain result
require 'smear'
--------------------------------------------------------------------



-- Subroutines------------------------------------------------------
function ape(alpha)
  local c = {}
  local ndims = #qopqdp.lattice()
  local norm = 0.5 * ndims * (ndims - 1)
  for mu = 1, ndims do
    c[mu] = { [0] = 1 - alpha }
    for nu = 1, ndims do
      if nu ~= mu then
        c[mu][nu] = alpha / norm
        c[mu][-nu] = alpha / norm
      end
    end
  end
  return c
end

-- Smearing only spatial links
function ape_space(alpha)
  local c = {}
  local ndims = #qopqdp.lattice()
  local norm = 0.5 * ndims * (ndims - 1) - 2
  for mu = 1, ndims - 1 do
    c[mu] = { [0] = 1 - alpha }
    for nu = 1, ndims - 1 do
      if nu ~= mu then
        c[mu][nu] = alpha / norm
        c[mu][-nu] = alpha / norm
      end
    end
  end
  return c
end

-- Return space--space and space--time plaquettes
-- Normalizations: (ndims - 1) * volume for space--time
--                 0.5 * (ndims - 2) * (ndims - 1) * volume for space--space
function getplaq(g)
  local ndims = #qopqdp.lattice()
  local norm_t = (ndims - 1) * g:lattice():volume()
  local norm_s = 0.5 * (ndims - 2) * norm_t
  local ss, st = g:action({plaq=1})
  ss = ss / norm_s
  st = st / norm_t
  return ss, st, 0.5 * (ss + st)
end
--------------------------------------------------------------------



-- Actual tests----------------------------------------------------
nx = 4
nt = 4
ls = { nx, nx, nx, nt }
L = qopqdp.lattice(ls)

TESTON()    -- Don't comment tmp-->out from here to TESTOFF

seed = 987654321
qopqdp.seed(seed)

-- Set up gauge field by some heatbath updates on a random configuration
-- Reunitarize at the end
local nrep = 10
local nhb = 1
local nor = 1
local beta = 6.0
local coeffs = {plaq=1}
g = qopqdp.gauge()
g:random()
g:heatbath(nrep, nhb, nor, beta, coeffs)
local devavg, devmax = g:checkSU()
printf("orig unitarity deviation avg: %.4g  max: %.4g\n", devavg, devmax)
g:makeSU()
devavg, devmax = g:checkSU()
printf("new  unitarity deviation avg: %.4g  max: %.4g\n", devavg, devmax)

-- Do various smearings and check chain rule
-- Also monitor plaquette as simple sanity check on smearing
-- First print unsmeared plaquette
local ss, st, tot = getplaq(g)
printf("Orig_plaq ss: %.8g  st: %.8g  tot: %.8g\n", ss, st, tot)

-- First smearing: 4d APE
smear = {}
local alpha = 0.5
smear[1] = { type="staples", coeffs = ape(alpha) }
myprint("smear = ", smear, "\n")
-- smearGauge expects an action object in the first argument...
local sg = smearGauge({g = g}, smear)

-- Check plaquette
--ss, st, tot = getplaq(sg)
printf("APE_plaq  ss: %.8g  st: %.8g  tot: %.8g\n", ss, st, tot)

-- Check derivative
-- TODO...

-- Second smearing: BSM-style HYP
alpha = {0.4, 0.5, 0.5}
smear[1] = { type="hyp", alpha = alpha }
sg = smearGauge({g = g}, smear)

-- Check plaquette
ss, st, tot = getplaq(sg)
printf("HYP_plaq  ss: %.8g  st: %.8g  tot: %.8g\n", ss, st, tot)

-- Check derivative
-- TODO...

-- Third smearing: adjoint rep
-- TODO...

TESTOFF()   -- Resume commenting tmp-->out
--------------------------------------------------------------------
