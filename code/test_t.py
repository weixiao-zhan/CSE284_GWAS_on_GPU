from scipy.stats import t
import numpy as np
import time


stat = 100*(np.random.random(10_000)-0.5)
df = 100
print(stat.shape)

start = time.time()
result_batch = t.sf(stat, df)
end = time.time()
print(end - start)

result_iter = np.zeros_like(result_batch)
start = time.time()
for i in range(stat.shape[0]):
    result_iter[i] = t.sf(stat[i], df)
end = time.time()
print(end - start)