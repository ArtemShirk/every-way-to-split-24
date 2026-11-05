"""Compute the first K nontrivial zeta zeros (imaginary parts) with mpmath. Usage: zeros_calc.py START STOP STEP OUT"""
import sys, mpmath, json
mpmath.mp.dps = 15
if len(sys.argv) < 5:
    print(__doc__)
    print('Example: python3 zeros_calc.py 1 2000 1 zeros_1_2000.json   (the zero files are already included)')
    sys.exit()
a, b, s, out = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
res = {}
for k in range(a, b + 1, s):
    res[k] = float(mpmath.zetazero(k).imag)
json.dump(res, open(out, 'w'))
