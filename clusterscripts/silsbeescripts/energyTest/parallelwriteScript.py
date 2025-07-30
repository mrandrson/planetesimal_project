import numpy as np
import math
import subprocess
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from multiprocessing import cpu_count

def run_job(args):
    subprocess.run([
        "python", "ShuOutput.py",
        str(args[0]), str(args[1]), str(args[2]), str(args[3]),
        str(args[4]), str(args[5]), str(args[6]), str(args[7])
    ], check=True)

if __name__ == "__main__":
    cs = math.sqrt(3.36) * 100
    #bvec = np.logspace(9, 17, 10)
    #Nv = 10
    #vvec = np.logspace(np.log10(cs / 10), np.log10(cs * 5), Nv)
    bvec = np.logspace(9, 18, 30)
    Nv = 10
    vvec = np.logspace(np.log10(cs / 100), np.log10(cs*5), Nv)
    jobs = [(bvec[i], bvec[i+1], vvec[j], 1.0829e13, i, j, 1, 100)
            for i in range(len(bvec)-1) for j in range(Nv)]

    workers = cpu_count()
    total = len(jobs)
    print(f"Launching {total} jobs on up to {workers} cores…")

    start = time.perf_counter()
    completed = 0
    with ProcessPoolExecutor(max_workers=workers) as exe:
        futures = {exe.submit(run_job, job): job for job in jobs}
        for fut in as_completed(futures):
            completed += 1
            elapsed = time.perf_counter() - start
            avg = elapsed / completed
            rem = avg * (total - completed)
            print(f"[{completed}/{total}] elapsed {elapsed:.1f}s, ETA {rem/60:.1f}m")

    total_time = time.perf_counter() - start
    print(f"All {total} jobs completed in {total_time/60:.2f} minutes ({total_time:.1f}s).")

