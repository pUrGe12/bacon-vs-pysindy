"""Convert ODEBench's equation list into data/odebench_equations.json (run once).

ODEBench (d'Ascoli et al., "ODEFormer", ICLR 2024) is 63 ODE systems, mostly from
Strogatz's textbook, each with fixed constants and two initial conditions.

    curl -sL -o /some/new/dir/strogatz_equations.py \
        https://raw.githubusercontent.com/sdascoli/odeformer/main/odeformer/odebench/strogatz_equations.py
    .venv/bin/python import_odebench.py /some/new/dir/strogatz_equations.py

The source file is read with `ast`, never imported or executed. Only literal values
(and constant arithmetic such as 8/3) are accepted.
"""
import ast
import hashlib
import json
import operator
import sys

SOURCE_URL = ("https://raw.githubusercontent.com/sdascoli/odeformer/main/"
              "odeformer/odebench/strogatz_equations.py")
# ODEBench's own trajectory settings (GPBench/ODEBench, generate_trajectories.py):
# t_span (0, 10), t_eval = np.linspace(0, 10, 150).
T_START, T_END, N_SAMPLES = 0.0, 10.0, 150

_OPS = {ast.Add: operator.add, ast.Sub: operator.sub, ast.Mult: operator.mul,
        ast.Div: operator.truediv, ast.USub: operator.neg, ast.UAdd: operator.pos}


def _value(node):
    if isinstance(node, ast.Constant):
        return node.value
    if isinstance(node, ast.List) or isinstance(node, ast.Tuple):
        return [_value(e) for e in node.elts]
    if isinstance(node, ast.Dict):
        return {_value(k): _value(v) for k, v in zip(node.keys, node.values)}
    if isinstance(node, ast.BinOp) and type(node.op) in _OPS:
        return _OPS[type(node.op)](_value(node.left), _value(node.right))
    if isinstance(node, ast.UnaryOp) and type(node.op) in _OPS:
        return _OPS[type(node.op)](_value(node.operand))
    raise ValueError(f"not a literal: {ast.dump(node)[:80]}")


def main(path):
    raw = open(path, "rb").read()
    tree = ast.parse(raw)
    eqs = None
    for stmt in tree.body:
        if isinstance(stmt, ast.Assign) and any(getattr(t, "id", "") == "equations"
                                                for t in stmt.targets):
            eqs = _value(stmt.value)
    if eqs is None:
        raise SystemExit("no `equations = [...]` found")
    out = {
        "source_url": SOURCE_URL,
        "source_sha256": hashlib.sha256(raw).hexdigest(),
        "time": {"t_start": T_START, "t_end": T_END, "n_samples": N_SAMPLES,
                 "from": "GPBench/ODEBench generate_trajectories.py config"},
        "equations": eqs,
    }
    with open("data/odebench_equations.json", "w") as fh:
        json.dump(out, fh, indent=1)
    dims = {}
    for e in eqs:
        dims[e["dim"]] = dims.get(e["dim"], 0) + 1
    print(f"{len(eqs)} systems written; by dimension: {dict(sorted(dims.items()))}")


if __name__ == "__main__":
    main(sys.argv[1])
