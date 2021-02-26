#!/bin/bash
## Project:
#SBATCH --account=NN9526K
##SBATCH --qos=preproc
## Job name:
#SBATCH --job-name=TGF-lin
## Wall time limit:
#SBATCH --time=0-12:00:0
## Number of nodes and task er node
#SBATCH --nodes=8 --ntasks-per-node=1 --cpus-per-task=32
#SBATCH --mail-user=david.sarria@uib.no

# asks SLURM to send the USR1 signal 120 seconds before end of the time limit
# in case we want to recover the data if the job is terminated due to time limit
#SBATCH --signal=B:USR1@120
# define the handler function
# note that this is not executed here, but when the associated signal is sent

ALT=9  
BEAM=15  
WKDR=/cluster/work/users/${USER}/

save_function_timeout()
{
    echo "function save_function_timeout called at $(date)"
    cat ${SCRATCH}/output_ascii/*.out > ${SCRATCH}/fused_${SLURM_JOB_ID}.out
    mkdir -p ${WKDR}/SIMULATION_DATAFILES/Geant4_TGF_propagation/all_data/${ALT}_${BEAM}
    cp ${SCRATCH}/fused_${SLURM_JOB_ID}.out ${WKDR}/SIMULATION_DATAFILES/Geant4_TGF_propagation/all_data/${ALT}_${BEAM}
}
# call save_function once we receive USR1 signal
trap 'save_function_timeout' USR1

save_function()
{
    echo "function save_function called at $(date)"
    cat ${SCRATCH}/output_ascii/*.out > ${SCRATCH}/fused_${SLURM_JOB_ID}.out
    mkdir -p ${WKDR}/SIMULATION_DATAFILES/Geant4_TGF_propagation/all_data/${ALT}_${BEAM}
    cp ${SCRATCH}/fused_${SLURM_JOB_ID}.out ${WKDR}/SIMULATION_DATAFILES/Geant4_TGF_propagation/all_data/${ALT}_${BEAM}
}
##

## Recommended safety settings:
set -o errexit # Make bash exit on any error
set -o nounset # Treat unset variables as errors

## Software modules

module restore system   # Restore loaded modules to the default
module load foss/2020b
module load icc/2019.1.144-GCC-8.2.0-2.31.1
module load CMake/3.9.1
module load Python/3.6.6-intel-2018b
module load iimpi/2020b 
module load Boost/1.68.0-intel-2018b-Python-3.6.6
module load iccifort/2020.4.304
module load impi/2019.9.304-iccifort-2020.4.304 
module load Python/3.6.6-intel-2018b
module load intel/2020a 


# module list             # List loaded modules, for easier debugging

## Run the application
echo "Started at $(date)"

cd ${SLURM_SUBMIT_DIR}
echo ${SLURM_NTASKS}
echo ${SCRATCH}
cd ${SCRATCH}

cp ${SLURM_SUBMIT_DIR}/TGF_Propa ${SCRATCH}
cp ${SLURM_SUBMIT_DIR}/vis.mac ${SCRATCH}
cp -R ${SLURM_SUBMIT_DIR}/TLE_data/ ${SCRATCH}
cp -R ${SLURM_SUBMIT_DIR}/mag_data/ ${SCRATCH}
mkdir output_ascii

# the "&" after the compute steps and "wait" are important
echo "Running Geant4"

a=1
b=3
rand_nb=$((a+RANDOM%(b-a))).$((RANDOM%999))

START=0
for (( ii=${START}; ii<${SLURM_JOB_NUM_NODES}; ii++ ))
do
   sleep $rand_nb
   srun -n1 -N1 ./TGF_Propa 800000000 ${ALT} 0.0 0.0 0 ${BEAM} 0 Gaussian 400 &
done
# n1 = 1 task
# N1 = 1 node
# "&" is important

# argument list:
# number_st = argv[1];
# settings.SOURCE_ALT = std::stod(argv[2]);
# settings.SOURCE_LAT = std::stod(argv[3]);
# settings.SOURCE_LONG = std::stod(argv[4]);
# settings.SOURCE_SIGMA_TIME = std::stod(argv[5]);
# settings.SOURCE_OPENING_ANGLE = std::stod(argv[6]);
# settings.TILT_ANGLE = std::stod(argv[7]);
# settings.BEAMING_TYPE = argv[8];
# settings.record_altitude = std::stod(argv[9]);

wait

echo "End Geant4 run"
save_function
wait

exit 0

