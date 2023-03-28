import sys
sys.path.insert(0, '../')
import modules as mod
import json
import time

method = 'pairwise'
w = 'w'
n_ary = 'RR'

# No trim
start = time.perf_counter()
lib = mod.SimilarityCalculator(trim_frac=None, n_ary=n_ary)
method_func = getattr(lib, f'calculate_{method}')
new_sims = method_func()

with open(f'{w}_{method}_{n_ary}.txt', 'w') as file:
    file.write(json.dumps(new_sims, indent=4))

end = time.perf_counter()
with open("benchmark.txt", "a") as f:
    f.write(f'{w}_{method}_{n_ary}: Finished in {round(end-start,2)} second\n\n')

# 10% outlier trim
trim_frac = 0.1
start = time.perf_counter()
lib = mod.SimilarityCalculator(trim_frac=trim_frac, n_ary=n_ary)
method_func = getattr(lib, f'calculate_{method}')
new_sims = method_func()

with open(f'{w}_{method}_{n_ary}_t{int(float(trim_frac) * 100)}.txt', 'w') as file:
    file.write(json.dumps(new_sims, indent=4))

end = time.perf_counter()
with open("benchmark.txt", "a") as f:
    f.write(f'{w}_{method}_{n_ary}_t{int(float(trim_frac) * 100)}: Finished in {round(end-start,2)} second\n\n')

# 20% outlier trim
trim_frac = 0.2
start = time.perf_counter()
lib = mod.SimilarityCalculator(trim_frac=trim_frac, n_ary=n_ary)
method_func = getattr(lib, f'calculate_{method}')
new_sims = method_func()

with open(f'{w}_{method}_{n_ary}_t{int(float(trim_frac) * 100)}.txt', 'w') as file:
    file.write(json.dumps(new_sims, indent=4))

end = time.perf_counter()
with open("benchmark.txt", "a") as f:
    f.write(f'{w}_{method}_{n_ary}_t{int(float(trim_frac) * 100)}: Finished in {round(end-start,2)} second\n\n')

# 30% outlier trim
trim_frac = 0.3
start = time.perf_counter()
lib = mod.SimilarityCalculator(trim_frac=trim_frac,n_ary=n_ary)
method_func = getattr(lib, f'calculate_{method}')
new_sims = method_func()

with open(f'{w}_{method}_{n_ary}_t{int(float(trim_frac) * 100)}.txt', 'w') as file:
    file.write(json.dumps(new_sims, indent=4))

end = time.perf_counter()
with open("benchmark.txt", "a") as f:
    f.write(f'{w}_{method}_{n_ary}_t{int(float(trim_frac) * 100)}: Finished in {round(end-start,2)} second\n\n')

