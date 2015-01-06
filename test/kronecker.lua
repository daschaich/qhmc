--------------------------------------------------------------------
-- This example is part of the regression test suite.  Functions starting
-- with TEST are part of the test framework and can be ignored.

-- Check simple 2x2 Kronecker product
--                     /  0  5  0 10 \
-- / 1 2 \ . / 0 5 \ = |  6  7 12 14 |
-- \ 3 4 /   \ 6 7 /   |  0 15  0 20 |
--                     \ 18 21 24 28 /

require 'gaugeObservables'  -- For plaquette
require 'smear'             -- "kronecker" is a type of smearing!

-- TODO: Check (2d?) adjoint plaquette against MILC
--------------------------------------------------------------------



--------------------------------------------------------------------
Nc = 2
nx = 4
nt = 4
ls = { nx, nx, nx, nt }
qopqdp.defaultNc(Nc)
L = qopqdp.lattice(ls)

TESTON()    -- Don't comment tmp-->out from here to TESTOFF

seed = 987654321
qopqdp.seed(seed)

-- Set up two gauge fields by some heatbath updates on random configurations
-- Reunitarize each and check plaquette
local nrep = 10
local nhb = 1
local nor = 1
local beta = 4.0
local coeffs = {plaq=1}
printf("Gauge field A:\n")
g1 = L:gauge()
g1:random()
g1:heatbath(nrep, nhb, nor, beta, coeffs)
printf("orig unitarity deviation avg: %.4g  max: %.4g\n", g1:checkSU())
g1:makeSU()
printf("new  unitarity deviation avg: %.4g  max: %.4g\n", g1:checkSU())
printf("fund plaq  ss: %.8g  st: %.8g  tot: %.8g\n", plaq(g1))

printf("Gauge field B:\n")
g2 = L:gauge()
g2:random()
g2:heatbath(nrep, nhb, nor, beta, coeffs)
printf("orig unitarity deviation avg: %.4g  max: %.4g\n", g2:checkSU())
g2:makeSU()
printf("new  unitarity deviation avg: %.4g  max: %.4g\n", g2:checkSU())
printf("fund plaq  ss: %.8g  st: %.8g  tot: %.8g\n", plaq(g2))

-- TODO: Do Kronecker product and check product plaquette here...
smear = {}
smear[1] = { type="kronecker", adj = {"false", "false"} }
local sg = smearGauge({g1 = g1, g2 = g2}, smear)
printf("prod plaq  ss: %.8g  st: %.8g  tot: %.8g\n", plaq(sg))


-- Set x-links at first sites to test case described above
-- Gauge field directions are 1, 2, 3, 4!
--[[
test = g2(1):point({0, 0, 0, 0}, 0, 0)
myprint("", test, "\n")
toSet = g2(1):point({0, 0, 0, 0}, 0, 1) -- Picks up complex pointer type...
g2(1):point({0, 0, 0, 0}, 0, 0, toSet)
test = g2(1):point({0, 0, 0, 0}, 0, 0)
myprint("", test, "\n")
]]
g1(1):point({0, 0, 0, 0}, 0, 0, 1)
g1(1):point({0, 0, 0, 0}, 0, 1, 2)
g1(1):point({0, 0, 0, 0}, 1, 0, 3)
g1(1):point({0, 0, 0, 0}, 1, 1, 4)
g2(1):point({0, 0, 0, 0}, 0, 0, 0)
g2(1):point({0, 0, 0, 0}, 0, 1, 5)
g2(1):point({0, 0, 0, 0}, 1, 0, 6)
g2(1):point({0, 0, 0, 0}, 1, 1, 7)

-- Do Kronecker product again and check resulting x-link at first site
-- TODO...
--local sg = smearGauge({g1 = g1, g2 = g2}, smear)


TESTOFF()   -- Resume commenting tmp-->out
--------------------------------------------------------------------
