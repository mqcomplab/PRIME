import argparse
import modules as mod

"""Add an argument for trim (optional) and Call the `gen_method_max` with the parsed trim argument
No trim
>>> python rep.py
Trim 10% outliers
>>> python rep.py -t 0.1
"""
parser = argparse.ArgumentParser(description='Generate method max with optional trim')
parser.add_argument('-t', '--trim', type=float, default=None,
                    help='Trim parameter for gen_method_max method')
args = parser.parse_args()


if args.trim is None:
    mod.gen_method_max()
else:
    mod.gen_method_max(trim=args.trim)