#!/bin/bash

# Parameters
SETUP=1

if [[ $SETUP == 1 ]]; then
  # Data distribution setting
  SETTING_LIST=(1)
  # Survival model types
  SURV_MODEL_TYPE_LIST=("grf") #  "cox"
  # Censoring model types
  CENS_MODEL_TYPE_LIST=("grf") #  "cox"
  # List of numbers of features
  N_FEAT_LIST=(20)
  # List of training sample sizes
  N_TRAIN_LIST=(200 500 1000 2000)
  # List of censoring training sample sizes
  N_TRAIN_CENS_LIST=(100)
  # List of calibration sample sizes
  N_CAL_LIST=(200)
  # Sequence of batches for parallel simulation
  BATCH_LIST=$(seq 1 1)

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
    for N_FEAT in "${N_FEAT_LIST[@]}"; do
      for N_CAL in "${N_CAL_LIST[@]}"; do
        for N_TRAIN in "${N_TRAIN_LIST[@]}"; do
          for N_TRAIN_CENS in "${N_TRAIN_CENS_LIST[@]}"; do
            for SURV_MODEL_TYPE in "${SURV_MODEL_TYPE_LIST[@]}"; do
              for CENS_MODEL_TYPE in "${CENS_MODEL_TYPE_LIST[@]}"; do

                # Generate a unique and interpretable file name based on the input parameters
                JOBN="setup_${SETUP}/setting${SETTING}_surv_${SURV_MODEL_TYPE}_cens_${CENS_MODEL_TYPE}_feat${N_FEAT}_train${N_TRAIN}_trainc${N_TRAIN_CENS}_cal${N_CAL}_batch${BATCH}.txt"
                OUT_FILE=$OUT_DIR"/"$JOBN
                #ls $OUT_FILE
                COMPLETE=0

                if [[ -f $OUT_FILE ]]; then
                  COMPLETE=1
                fi

                if [[ $COMPLETE -eq 0 ]]; then
                  # R script to be run with command line arguments
                  SCRIPT="./experiment_1.sh $SETUP $SETTING $SURV_MODEL_TYPE $CENS_MODEL_TYPE $N_FEAT $N_TRAIN $N_TRAIN_CENS $N_CAL $BATCH"

                  # Define job name for this configuration
                  OUTF=$LOGS"/"$JOBN".out"
                  ERRF=$LOGS"/"$JOBN".err"

                  # Assemble slurm order for this job
                  ORD=$ORDP" -J "$JOBN" -o "$OUTF" -e "$ERRF" $SCRIPT"

                  # Print order
                  echo $ORD
                  # Submit order
                  $ORD
                  # Run command now
                  #./$SCRIPT
                  
                fi

              done
            done
          done
        done
      done
    done
  done
done
