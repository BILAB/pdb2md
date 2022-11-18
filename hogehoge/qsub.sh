#!/bin/zsh
DIR=(
	D5SLU6
)
for i in $DIR
do
    echo $1
    cd $i/amber
    / usr/local/bin/qsub ./totalrun.sh - N $i
    cd ../.. /
done