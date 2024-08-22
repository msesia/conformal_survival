#!/bin/bash

# Parameters
SETUP=1

if [[ $SETUP == 1 ]]; then
  # Data distribution setting
  SETTING_LIST=(3)
  # List of calibration sample sizes
  N_CAL_LIST=(500)
  # List of training sample sizes (interpreted as num_samples_train)
  N_TRAIN_LIST=(1000)
  # Sequence of batches for parallel simulation
  BATCH_LIST=$(seq 1 10)

  MEMO=5G
fi

# Slurm parameters
TIME=00-00:20:00                    # Time required (20 m)
CORE=1                              # Cores required (1)

# Assemble order prefix
ORDP="sbatch --mem="$MEMO" --nodes=1 --ntasks=1 --cpus-per-task=1 --time="$TIME" --account=sesia_1124 --partition=main"

# Create directory for log files
LOGS="logs"
mkdir -p $LOGS
mkdir -p $LOGS"/setup_"$SETUP

OUT_DIR="results"
mkdir -p $OUT_DIR
mkdir -p $OUT_DIR"/setup_"$SETUP

# Loop over configurations
for BATCH in $BATCH_LIST; do
  for SETTING in "${SETTING_LIST[@]}"; do
    for N_CAL in "${N_CAL_LIST[@]}"; do
      for N_TRAIN in "${N_TRAIN_LIST[@]}"; do

        # Generate a unique and interpretable file name based on the input parameters
        JOBN="setup_${SETUP}/setting${SETTING}_train${N_TRAIN}_cal${N_CAL}_batch${BATCH}.txt"
        OUT_FILE=$OUT_DIR"/"$JOBN
        #ls $OUT_FILE
        COMPLETE=0

        if [[ -f $OUT_FILE ]]; then
          COMPLETE=1
        fi

        if [[ $COMPLETE -eq 0 ]]; then
          # R script to be run with command line arguments
          SCRIPT="./experiment_1.sh $SETUP $SETTING $N_TRAIN $N_CAL $BATCH"

          # Define job name for this configuration
          OUTF=$LOGS"/"$JOBN".out"
          ERRF=$LOGS"/"$JOBN".err"

          # Assemble slurm order for this job
          ORD=$ORDP" -J "$JOBN" -o "$OUTF" -e "$ERRF" $SCRIPT"

          # Print order
          echo $ORD
          # Submit order
          #$ORD
          # Run command now
          #./$SCRIPT

          
        fi

      done
    done
  done
done