import numpy as np
import argparse
parser = argparse.ArgumentParser(description='Create a jobs file for build_stamps')
parser.add_argument("--max_idx", action  ="store",
                                 dest    = "max_idx",
                                 default = 1000) #total for DR is 212404
parser.add_argument("--num_proc", action = "store",
                                  dest   = "num_proc",
                                  default = 10)
parser.add_argument('--version', action='version', version='0.1')
results = parser.parse_args()

if __name__=="__main__":
    idx    = np.arange(int(results.max_idx), dtype=int)
    idxmat = idx.reshape((int(results.num_proc), -1))
    np.savetxt("jobs.txt", idxmat, fmt="%d ")

