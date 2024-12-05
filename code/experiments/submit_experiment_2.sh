#!/bin/bash

# Parameters
SETUP=0

if [[ $SETUP == 0 ]]; then
  # Data distribution setting
  DATA_LIST=("VALCT" "PBC" "GBSG" "METABRIC")
  # Survival model types
  SURV_MODEL_TYPE_LIST=("grf" "rf" "cox" "survreg")
  #SURV_MODEL_TYPE_LIST=("cox")
  # Censoring model types
  CENS_MODEL_TYPE_LIST=("grf")
  # Subsampling for training set
  TRAIN_PROP_LIST=(1.0)
  # Sequence of batches for parallel simulation
  BATCH_LIST=$(seq 1 4)

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
mkdir -p $LOGS"/data"

OUT_DIR="results"
mkdir -p $OUT_DIR
mkdir -p $OUT_DIR"/data"

# Loop over configurations
for BATCH in $BATCH_LIST; do
  for DATA in "${DATA_LIST[@]}"; do
    for SURV_MODEL_TYPE in "${SURV_MODEL_TYPE_LIST[@]}"; do
      for CENS_MODEL_TYPE in "${CENS_MODEL_TYPE_LIST[@]}"; do
        for TRAIN_PROP in "${TRAIN_PROP_LIST[@]}"; do

          # Generate a unique and interpretable file name based on the input parameters
          JOBN="data/${DATA}_surv_${SURV_MODEL_TYPE}_cens_${CENS_MODEL_TYPE}_train_${TRAIN_PROP}_batch${BATCH}.txt"
          OUT_FILE=$OUT_DIR"/"$JOBN
          #ls $OUT_FILE
          COMPLETE=0

          if [[ -f $OUT_FILE ]]; then
            COMPLETE=1
          fi

          if [[ $COMPLETE -eq 0 ]]; then
            # R script to be run with command line arguments
            SCRIPT="./experiment_2.sh $DATA $SURV_MODEL_TYPE $CENS_MODEL_TYPE $TRAIN_PROP $BATCH"

            # Define job name for this configuration
            OUTF=$LOGS"/"$JOBN".out"
            ERRF=$LOGS"/"$JOBN".err"

            # Assemble slurm order for this job
            ORD=$ORDP" -J "$JOBN" -o "$OUTF" -e "$ERRF" $SCRIPT"

            # Print order
            echo $ORD
            # Submit order
#            $ORD
            # Run command now
            #./$SCRIPT
            
          fi

        done
      done
    done
  done
done
