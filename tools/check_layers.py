#!/usr/bin/env python3
"""Check the layer rule: no module outside src/reference/ may use a module inside it.  Exit 1 on a violation.
Prints the modules and the lines the program (app/main.f90) reaches."""
import glob, re, sys

def uses(f):
    return set(m.lower() for m in re.findall(r'^\s*use\s*(?:,\s*intrinsic\s*::)?\s*(\w+)', open(f, errors='ignore').read(), re.M | re.I))

mods = {}
for f in glob.glob('src/**/*.[fF]90', recursive=True):
    m = re.search(r'^\s*module\s+(\w+)\s*$', open(f, errors='ignore').read(), re.M | re.I)
    if m:
        mods[m.group(1).lower()] = f
ref = {m for m, f in mods.items() if '/reference/' in f}
bad = [(m, u) for m, f in mods.items() if m not in ref for u in uses(f) if u in ref]
for m, u in bad:
    print(f'violation: production module {m} uses reference module {u}')
seen, stack = set(), list(uses('app/main.f90'))
while stack:
    m = stack.pop()
    if m in seen or m not in mods:
        continue
    seen.add(m)
    stack += list(uses(mods[m]))
nl = sum(sum(1 for _ in open(mods[m], errors='ignore')) for m in seen)
print(f'program reaches {len(seen)} modules, {nl} lines; reference layer: {len(ref)} modules')
sys.exit(1 if bad or (seen & ref) else 0)
