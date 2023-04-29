import argparse
import modules as mod

"""Call the `gen_method_max` with the parsed trim and index argument
No trim
>>> python rep.py -i SM
Trim 10% outliers
>>> python rep.py -t 0.1 -i RR
"""
parser = argparse.ArgumentParser(description='Generate method max with optional trim and n_ary')
parser.add_argument('-t', '--trim', type=float, default=None,
                    help='Trim parameter for gen_method_max method')
parser.add_argument('-i', '--index', type=str, default='RR',
                    help='n_ary parameter for gen_method_max method')
args = parser.parse_args()

if args.trim is None:
    mod.gen_method_max(n_ary=args.index)
else:
    mod.gen_method_max(trim=args.trim, n_ary=args.index)