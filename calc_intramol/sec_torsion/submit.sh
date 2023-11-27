sbatch -p single --export=ALL,OMP_NUM_THREADS=8 -J B3LYP -N 1 -c 8 -t 54:10:00 --mem=2000 ./scan_start.sh
