"""Finite matrix checks of the leaking rotating-dephasing example.
The adaptive upper bound is proved by the common-processor argument in the text,
not by these finite tests.
"""
import numpy as np

def trace_distance(a, b):
    return np.abs(np.linalg.eigvalsh(a-b)).sum()/2

def tensor_power(a, n):
    result = np.array([[1.0]])
    for _ in range(n): result = np.kron(result, a)
    return result

zero = np.diag([1.0, 0.0])
for theta in (0.0, 0.01, 0.2, 0.7, np.pi/2):
    c, s = np.cos(theta), np.sin(theta)
    v = np.array([c, s]); w = np.array([-s, c])
    p, q = np.outer(v,v), np.outer(w,w)
    ideal = np.diag([0.5, 0, 0, 0.5])
    actual = (np.kron(p, zero) + np.kron(q, np.eye(2)-zero))/2
    assert np.isclose(trace_distance(actual, ideal), s)
    assert np.isclose(np.trace(np.diag([0,1,1,0]) @ actual), s*s)
    for n in range(1, 7):
        assert np.isclose(trace_distance(tensor_power(p,n), tensor_power(zero,n)),
                          np.sqrt(max(0, 1-c**(2*n))), atol=2e-8)
print("Rotating-dephasing example: Choi distance, leakage and n=1..6 product distances verified.")
