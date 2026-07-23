"""
singular_to_msolve.py — convert Singular-dumped generator files (one poly per
line, short=0 output: explicit '*' and '^') into msolve input format:

  line 1: comma-separated variable names
  line 2: field characteristic
  then the polynomials, comma-separated (one per line, trailing comma on all
  but the last).

Usage: python singular_to_msolve.py q1_gens.txt q1.ms
Variables are the fixed 31-variable set of the T5 dehomogenized probes.
"""
import sys

AV = [f'a{i}' for i in range(5)]
BV = [f'b{i}' for i in range(7)]
CV = [f'c{i}' for i in range(9)]
EV = [f'e{i}' for i in range(1, 11)]
ALLV = AV + BV + CV + EV
CHAR = 32003


def convert(src, dst):
    polys = []
    for line in open(src):
        p = line.strip()
        if not p or p == '0':
            continue
        polys.append(p)
    with open(dst, 'w') as f:
        f.write(','.join(ALLV) + '\n')
        f.write(f'{CHAR}\n')
        for i, p in enumerate(polys):
            f.write(p + (',\n' if i < len(polys) - 1 else '\n'))
    print(f'{dst}: {len(polys)} polynomials, {len(ALLV)} vars, char {CHAR}')


if __name__ == '__main__':
    convert(sys.argv[1], sys.argv[2])
