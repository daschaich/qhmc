require 'Lattice'
require 'Action'
require 'Evolver'
require 'topo'

--profile(1)
--verbosity(1)
--trace(1)

nx = 8
nt = 8
beta = 20

--L = Lattice{nx,nx,nt}
L = Lattice{nx,nx,nx,nt}
L:Seed(987654321)

G = L:GaugeField{nc=3}
G:Set("unit")

GA = Action{kind="gauge",style="plaquette",beta=beta,field=G}

E = Evolver{kind="heatbath",actions={GA},
	    nRepetitions=10,nHeatbath=1,nOverrelax=1}

function ploop()
  local nd = #L
  local pl = {}
  local i = nd-1
  local plpath = {}
  --for j=1,L[i] do plpath[j] = -i end -- changed convention
  for j=1,L[i] do plpath[j] = i end
  pl = G.field:loop(plpath)
  return pl
end

function topo()
  local tr,ti = pathDo(G.field, #L.latsize, paths0, coeffs0)
  return tr
end

local ps, pt, ptot = plaq()
printf("initial plaq: %g\n", ptot)
printf("initial ploop: %s\n", tostring(ploop()))

for i=1,10 do
  E:Run()
  ps, pt, ptot = plaq()
  printf("plaq: %g\n", ptot)
  printf("ploop: %s\n", tostring(ploop()))
  --printf("topo: %s\n", tostring(topo()))
end
