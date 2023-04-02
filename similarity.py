import argparse
import modules as mod
import json
import time

# Parse command-line arguments
parser = argparse.ArgumentParser()
parser.add_argument('-m', '--method', help='Method to use for similarity calculation. \
                    (pairwise, union, medoid, outlier)', required=True)
parser.add_argument('-n', '--n_clusters', type=int, help='Number of clusters for analysis', 
                    required=True)
parser.add_argument('-i', '--index', help='Similarity Index to use (e.g. RR or SM)', 
                    required=True)
parser.add_argument('-t', '--trim_frac', type=float, help='Fraction of outliers to trim. \
                    (e.g. 0.1, default: None)', default=None)
parser.add_argument('-w', '--weighted', help='Weighing clusters by frames it contains. \
                    (default: True)', default=True)
parser.add_argument('-d', '--cluster_folder', help='Location of the cluster files directory', 
                    default="new_clusters/")
parser.add_argument('-s', '--summary_file', help='Location of CPPTRAJ cluster summary file', 
                    default="../Cpptraj_linkage_sieve_eps_2/summary")
args = parser.parse_args()

# Calculate similarities
start = time.perf_counter()
lib = mod.SimilarityCalculator(cluster_folder=args.cluster_folder, summary_file=args.summary_file, 
                               n_clusters=args.n_clusters, trim_frac=args.trim_frac, 
                               n_ary=args.index, weighted=args.weighted)
method_func = getattr(lib, f'calculate_{args.method}')
new_sims = method_func()
if args.weighted is True:
    w = "w"
else:
    w = "nw"
if args.trim_frac:
    with open(f'{w}_{args.method}_{args.index}_t{int(float(args.trim_frac) * 100)}.txt', 
              'w') as file:
        file.write(json.dumps(new_sims, indent=4))
    end = time.perf_counter()
    print(f"{w}_{args.method}_{args.index}_t{int(float(args.trim_frac) * 100)}: \
          Finished in {round(end-start,2)} second")
else:
    with open(f'{w}_{args.method}_{args.index}.txt', 'w') as file:
        file.write(json.dumps(new_sims, indent=4))
    end = time.perf_counter()
    print(f"{w}_{args.method}_{args.index}: Finished in {round(end-start,2)} second")
    
# python similarity.py -m medoid -n 11 -i RR
