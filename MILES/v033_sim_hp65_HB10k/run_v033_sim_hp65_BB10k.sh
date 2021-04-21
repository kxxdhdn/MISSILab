#!/bin/bash -l
#SBATCH -J v033_sim_hp65_HB10k
#SBATCH -o /dsm/herschel1/nuages/dhu/PAHPedia/v033_sim_hp65_HB10k/v033_sim_hp65_BB10k.out
#SBATCH -e /dsm/herschel1/nuages/dhu/PAHPedia/v033_sim_hp65_HB10k/v033_sim_hp65_BB10k.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=100-01:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dangning.hu@outlook.com

echo "Starting the job on" `date`
echo "The shell is "$SHELL
echo "We are in the directory "$PWD
echo "We are on "$HOSTNAME

cd /dsm/herschel1/nuages/dhu/MISSILE/MILES/programs
srun -n1 /dsm/herschel1/nuages/dhu/MISSILE/MILES/programs/fitpar_HB

echo "Finishing the job on" `date`
exit
