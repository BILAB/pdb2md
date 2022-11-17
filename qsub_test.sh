#ディレクトリのリスト/a/nrjevと/b/frekと/c/fcrneaolを変数DIRに格納
DIR="/a/nrjev /b/frek /c/fcrneaol"
#変数DIRの中のディレクトリをループ処理
for i in $DIR
do
cd $i/amber/
qsub totallrun.sh -N vekm
cd ..
done
#次のディレクトリに移動