#!/bin/sh
#PBS -q default
#PBS -l nodes=1:ppn=16:gpus=1
#PBS -l walltime=72:00:00

test $PBS_O_WORKDIR && cd $PBS_O_WORKDIR
# run the environment module
if test -f /home/apps/Modules/init/profile.sh; then
    . /home/apps/Modules/init/profile.sh
    module load amber22
elif test -f /usr/local/Modules/init/profile.sh; then
    . /usr/local/Modules/init/profile.sh
    module load amber22
fi

### Write your qsub script from here.
echo `hostname`

# トポロジーファイルの指定
topfile="../../../top/leap.parm7"
# 再開させたいrst7ファイルを指定
rstfile="../001/md.rst7"

pmemd.cuda_SPFP.MPI -O \
    -i md.in \
    -o md.out \
    -p ${topfile} \
    -c ${rstfile} \
    -r md.rst7 \
    -ref ${topfile}

