import argparse
import ast
import sys
import numpy as np
from model import SIR, SIRX

MODEL_SIR = 0
MODEL_SIRX = 1

FIRST_ROW = [
        ["Seconds", "S", "I", "R"],
        ["Seconds", "S", "I", "R", "X"]]

def run_cli():
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--model', default="sir", choices=["sir", "sirx"])
    parser.add_argument('-p', '--parameters', required=True)
    parser.add_argument('-i', '--initial-conds', required=True)
    parser.add_argument('-t', '--time', type=int)

    args = parser.parse_args()
    p = np.array(ast.literal_eval(args.parameters))
    w = np.array(ast.literal_eval(args.initial_conds))
    pop = np.sum(w)
    w0 = w/pop
    model = MODEL_SIR if args.model == "sir" else MODEL_SIRX
    kwargs = {}

    if (args.time):
        kwargs["tf_days"] = args.time

    if (model == MODEL_SIR):
        cls = SIR
    else:
        cls = SIRX

    # Call the desired solver
    out = cls(p, w0).solve(**kwargs)

    # Multiply by the population
    out[:,1:] *= pop

    # TODO: Improve this
    print(",".join(FIRST_ROW[model]))
    np.savetxt(sys.stdout, out, delimiter=",")

if __name__ == "__main__":
    run_cli()
