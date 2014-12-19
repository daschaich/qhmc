--------------------------------------------------------------------
-- This example is part of the regression test suite.  Functions starting
-- with TEST are part of the test framework and can be ignored.

-- To test chain rule for smeared fields, we start with a gauge field X
-- and some chain field, calculate
--   Re Tr[ f(X) C ]
-- then we add some small value to one element of X, calculate
-- the finite difference and compare to the smearChain result
require 'gaugeObservables'
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
printf("orig unitarity deviation avg: %.4g  max: %.4g\n", g:checkSU())
g:makeSU()
printf("new  unitarity deviation avg: %.4g  max: %.4g\n", g:checkSU())

-- Do various smearings and check chain rule
-- Also monitor plaquette as simple sanity check on smearing
-- First print unsmeared plaquette
printf("Orig_plaq  ss: %.8g  st: %.8g  tot: %.8g\n", plaq(g))

-- First smearing: 4d APE
smear = {}
local alpha = 0.5
smear[1] = { type="staples", coeffs = ape(alpha) }
myprint("smear = ", smear, "\n")
-- smearGauge expects an action object in the first argument...
local sg = smearGauge({g = g}, smear)

-- Check plaquette
printf("APE_plaq   ss: %.8g  st: %.8g  tot: %.8g\n", plaq(sg))

-- Check derivative
-- TODO...

-- Second smearing: stout
alpha = 0.1
smear[1] = { type="stout", rho = alpha }
myprint("smear = ", smear, "\n")
sg = smearGauge({g = g}, smear)

-- Check plaquette
printf("stout_plaq ss: %.8g  st: %.8g  tot: %.8g\n", plaq(sg))

-- Check derivative
-- TODO...

-- Third smearing: BSM-style HYP
alpha = {0.4, 0.5, 0.5}
smear[1] = { type="hyp", alpha = alpha }
myprint("smear = ", smear, "\n")
sg = smearGauge({g = g}, smear)

-- Check plaquette
printf("HYP_plaq   ss: %.8g  st: %.8g  tot: %.8g\n", plaq(sg))

-- Check derivative
-- TODO...

-- Fourth smearing: adjoint rep
-- TODO...

TESTOFF()   -- Resume commenting tmp-->out
--------------------------------------------------------------------
