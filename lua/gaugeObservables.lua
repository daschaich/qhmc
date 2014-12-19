--------------------------------------------------------------------
-- Return space--space plaquettes, space--time plaquettes and their average
-- Normalizations: (ndims - 1) * volume for space--time
--                 0.5 * (ndims - 2) * (ndims - 1) * volume for space--space
function plaq(g)
  local ndims = #(g:lattice())
  local norm_t = (ndims - 1) * g:lattice():volume()
  local norm_s = 0.5 * (ndims - 2) * norm_t
  local ss, st = g:action({plaq=1})
  ss = ss / norm_s
  st = st / norm_t
  return ss, st, 0.5 * (ss + st)
end

function plaqE(g)
  local ss, st, tot = plaq(g)
  local nd = #(g:lattice())
  local e = 0.5 * nd * (nd - 1) * (2 - ss - st)
  return e
end
--------------------------------------------------------------------
