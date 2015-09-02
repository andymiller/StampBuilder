# StampBuilder
Builds small image patches of SDSS data.  Organizes corresponding spectra.

Requires: tractor, astrometry

This is the roadmap for building a dataset of stamps on Edison from a catalog
fits file - this is IO bound, so we'll set up everything the serial queue.

1. 
  ```
  export QDO_BATCH_PROFILE=edison-serial
  ```

2. Create a `jobs.txt` file.  This is a textfile of arguments - each line 
   corresponds to one machine executing the job `build_stamps.sh 0 1 2 3 ...`

   Change `build_jobs_file.py` to reflect the number of stamps you want to create
   and how many stamps per core, then run it
  ```
  python build_jobs_file.py --max_idx=1000 --num_proc=10
  ```
  which will create a jobs.txt file that will build 1000 stamps using 10 processors.

3. Set up a QDO task, `stamps`
  ```
  qdo load stamps jobs.txt
  ```

2. Execute the task stamps on N cores
  ```
  qdo launch stamps 10 --batchopts "-A cosmo" --walltime 1:00:00 --batchqueue serial --script ./build_stamps.sh
  ```
  notes:
  * This calls for 10 serial cores, using account "cosmo"
  * walltime is H:MM:SS
  * `build_stamps.sh 0 1 2 3` calls `build_dataset.py` on each of the indices
    (rows in jobs.txt file)
