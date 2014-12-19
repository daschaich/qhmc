-- This example is part of the regression test suite.  Functions starting
-- with TEST are part of the test framework and can be ignored.

require 'gaugeObservables'
require 'topo'

nx = 4
nt = 4
ls = { nx, nx, nx, nt }
L = qopqdp.lattice(ls)

TESTON()

seed = 987654321
qopqdp.seed(seed)

g = qopqdp.gauge()
g:unit()
printf("plaq ss: %-8g  st: %-8g  tot: %-8g\n", plaq(g))

g:random()
printf("plaq ss: %-8g  st: %-8g  tot: %-8g\n", plaq(g))

local nrep = 10
local nhb = 1
local nor = 1
local beta = 10
local coeffs = {plaq=1}
g:heatbath(nrep, nhb, nor, beta, coeffs)
printf("plaq ss: %-8g  st: %-8g  tot: %-8g\n", plaq(g))

t0 = clock()
se,sq = symmEQ(g)
t0 = clock() - t0
printf("#time: %.8g\n", t0)
printf("se: %.8g\n", se)
printf("sq: %.8g\n", sq)

t0 = clock()
se,sq = symmEQ(g,1)
t0 = clock() - t0
printf("#time: %.8g\n", t0)
printf("se: %.8g\n", se)
printf("sq: %.8g\n", sq)

t0 = clock()
se,sq = symmEQ(g,1,"timeslices")
t0 = clock() - t0
printf("#time: %.8g\n", t0)
for i=1,#se do
  printf("%i  se: %.7g  sq: %.7g\n", i, se[i], sq[i])
end

TESTOFF()
