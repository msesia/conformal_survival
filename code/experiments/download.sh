DATA=$1

mkdir -p results_hpc

rsync -auv sesia@discovery.usc.edu:/home1/sesia/Workspace/conformal_survival/code/experiments/results/* results_hpc/
