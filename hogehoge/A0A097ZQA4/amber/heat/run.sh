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

topfile="../../top/leap.parm7"
prev_rst7="../minimize/min2.rst7"

for (( i=1;i<10;i++ ));do
    echo "Start the ${i} cycle."
    #まだmd${i}ステップが終わっていなければ(md${i}.rst7が存在していないなら)実行
    if [ ! -f md${i}.rst7 ];then
        pmemd.cuda_SPFP.MPI -O \
            -i md${i}.in \
            -o md${i}.out \
            -p ${topfile} \
            -c ${prev_rst7} \
            -ref ${prev_rst7} \
            -r md${i}.rst7 \
            -x md${i}.nc \
            -inf md${i}.info || exit $?
    #すでにmd${i}ステップが終わっていれば次のサイクルへ移動
    else
        echo "md${i}.rst7 exists. Starts next cycle."
    fi
    prev_rst7="md${i}.rst7"
done
# 不要であればmd1~md9.ncを消去
rm md[1-9].nc
# ambpdbでmd9.rst7をmd9.pdbファイルに変換
ambpdb -p ${topfile} -c md9.rst7 > md9.pdb
