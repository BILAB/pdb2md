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
test ${PBS_NP} || PBS_NP=8
DO_PARALLEL="mpirun -np ${PBS_NP} --mca orte_base_help_aggregate 0"
topfile="../../top/leap.parm7"
rstfile="../../top/leap.rst7"

${DO_PARALLEL} pmemd.MPI -O \
    -i min1.in \
    -o min1.out \
    -p ${topfile} \
    -c ${rstfile} \
    -r min1.rst7 \
    -inf min1.info || exit $?
${DO_PARALLEL} pmemd.MPI -O \
    -i min2.in \
    -o min2.out \
    -p ${topfile} \
    -c min1.rst7 \
    -r min2.rst7 \
    -inf min2.info || exit $?

# ambpdbコマンドでmin2操作後の座標ファイルをpdb形式のファイルに変換する
ambpdb -p ${topfile} -c min2.rst7 > min2.pdb
