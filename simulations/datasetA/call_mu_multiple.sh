#!/bin/bash
#SBATCH --job-name=gillespie_master
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH -o main_job.out
#SBATCH -e main_job.err

#01. N
#02. M
#03. niter
#04. seed
#05. init => 1: zeros, 2: ordered, 3: shuffled
#06. kT
#07. spring_k
#08. spring_l (-1: sequential spring_l)
#09. offset
#10. Vmax (-1: original Vmax)
#11. Drive (kT)
#12. Tau free (s) DNA free
#13. Tau bound (s) DNA bound


#14. file
#15. complete output ('y' = yes, 'n' = no)
#16. output trajectory ('y' = yes, 'n' = no)
#17. output time ('y' = yes, 'n' = no)
#18. output cycle ('y' = yes, 'n' = no)
#19. number of repetitions (only for complete = 'n')


c=0

mu=50

for k in {1..10}
do
    cp base_script.sh script_"$c".sh
    echo srun -p LP python3 10_gillespie_trajectories_script.py 30 6 5000000 "$c" 2 1 "$k"00 -1 0 -1 "$mu"0 1 1 mu"$mu"0_k"$k"00 n y y y 10 >> script_"$c".sh
    sbatch -p LP script_"$c".sh
    c=$(($c + 1))
done
