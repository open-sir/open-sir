# pylint: disable=C0103
""" Open-SIR CLI implementation """
import argparse
import ast
import sys
import numpy as np
from model import SIR, SIRX

MODEL_SIR = 0
MODEL_SIRX = 1

def run_cli():
    """ This function runs the main CLI routine """
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

    if args.time:
        kwargs["tf_days"] = args.time

    if model == MODEL_SIR:
        sol = SIR(p, w0, pop)
    else:
        sol = SIRX(p, w0, pop)

    sol.solve(**kwargs).export(sys.stdout)

if __name__ == "__main__":
    run_cli()
