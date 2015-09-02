# StampBuilder
Builds small image patches of SDSS data.  Organizes corresponding spectra.

Requires: tractor, astrometry


## ROADMAP
The following is an overview of the steps to build a dataset of stamps on Edison from a catalog
fits file - this is IO bound, so we'll set up everything the serial queue.  These scripts 
are set up to take a specific fits source file and create stamps output into
a fixed directory (specified in `build_stamps.sh`). 

  * source file:  `PhotoSpecBoss_andrewcmiller.fit`
  * output directory: `/project/projectdirs/das/acmiller/stamps`

Firstly - an example [casjob](http://skyserver.sdss.org/CasJobs/) query that this works with is (in context DR9)
```
select * into mydb.PhotoSpecBoss from specphotoall 
where 
  SciencePrimary=1 and
  (ra between 0 and 180) and 
  (dec between -18 and 18) and 
  survey = 'boss'
```
which can be downloaded to create a file `PhotoSpecBoss_andrewcmiller.fit`.  

Then, run the pipeline on NERSC with the following steps: 


1. First, edison-serial queue is special for QDO.
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
  which will create a `jobs.txt` file that will build 1000 stamps using 10 processors.

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
  * walltime is `H:MM:SS`
  * script `build_stamps.sh 0 1 2 3` calls `build_stamp.py` on each of the indices
    (rows in jobs.txt file)


