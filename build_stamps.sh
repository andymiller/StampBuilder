#!/bin/bash
for i in $*; do
  echo "STAMP IDX: $i"
  python build_stamp.py \
    --source_file=fits_files/PhotoSpecBoss_andrewcmiller.fit \
    --output_dir=/project/projectdirs/das/acmiller/stamps \
    --num_proc=1 \
    --idx_keep=$i
done
