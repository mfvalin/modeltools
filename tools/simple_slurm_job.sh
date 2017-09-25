#!/home/valin/dummy_shell_slurm
#SBATCH --time=00:01:00
#SBATCH --account=def-example
#SBATCH --output=/home/valin/listings/graham/%j.o
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=4
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
. /etc/profile
set -x
echo 'Hello, world!'
echo ================================================================
env | sort
echo ================================================================
ps -ef | grep -v root
echo ================================================================
lsb_release -a
echo ================================================================
which ifort
ifort --version
which mpif90
mpif90 --version
echo ================================================================
mpirun -h
echo ================================================================
mpirun -n 8 --npernode 4 /bin/hostname
sleep 30

