sbatch -p single --export=ALL,OMP_NUM_THREADS=8 -J mp2 -N 1 -c 8 -t 70:00:00 --mem=2000 ./scan_start.sh
