"""Exact finite bounds for the displayed smoothed Choi-model example.

Uses the manuscript's analytic second-order/Choi inequalities and pi < 22/7.
This does not reconstruct the finite model or prove those analytic inequalities.
"""
from fractions import Fraction as F
import math

k, M = 16, 2120
m = (M - 1) * (M - 2) // 2 + 1
# Rational upper bounds for all square roots in Upsilon_k.
roots = [(F(21, 2), F(3241, 1000)), (F(21, 4), F(2292, 1000)),
         (F(11), F(3317, 1000)), (F(2), F(1415, 1000)),
         (F(6), F(245, 100)), (F(3), F(1733, 1000))]
assert all(bound > 0 and bound**2 > radicand for radicand, bound in roots)
a, b, c, d, e, f = [bound for _, bound in roots]
upsilon = (a + 2*b + c + F(1, 4**k)*(2 + 2*d + e + 2*f)) / 104
phase = F(22, 7)**2 * (F(k, 8*M*m) + F(k*k, 2*m*m))
error = F(41, 32 * 4**k) + upsilon * phase
assert error < F(1, 10**9)
# Exact rounding: 22 < log2(M^2) < 23; no floating-point ceil is needed.
assert 2**22 < M*M < 2**23
assert k*(4*M+4) + 22 == 135766
assert k*(4*M+4) + 23 == 135767
print(f"Choi upper bound < {float(error):.12g}; constructed log2 dimension "
      f"{k*(4*M+4)+2*math.log2(M):.9f}; exact qubit ceiling 135767.")

# Certified margins for the 135 M^(log2 107) area bound. In each inequality
# the ratio decreases with p, so p >= 337/50 suffices; all comparisons below
# are integers after raising to the fiftieth power.
assert 107**50 > 2**337
th_squared_base = (750*135*68)**2 * 624
assert th_squared_base**50 * 9**337 * 16**724 < 625**724  # c_H > .0016
kh_squared_numerator = (10**8*135*2**17)**2 * 16
assert kh_squared_numerator**50 * 9**337 * 17**724 < 15**50 * 6250**724  # c_a > .00017
p=math.log2(107)
beta=1/(2*p+1)
th=750*135*3**p*68*math.sqrt(624)
kh=10**8*135*3**p*2**17
ca=(16*kh**(2*beta))**-1*(16/15)**(-beta)
assert 2.831e11 < th < 2.833e11 and 2.07e41 < (8*ca)**(-1/beta) < 2.08e41
print("Exact 107-profile prefactors c_H > .0016, c_a > .00017 passed; rounded displays checked.")
