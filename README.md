# pdb2md
## Overview
・指定したID(uniplot / PDB)の結晶構造(uniplotの場合はAlphaFoldDBの予測構造)の取得

・(for PDB)ミッシング領域の構築

・ディスオーダー領域の削除

・テンプレート構造を基にしたヘテロ原子や基質の挿入

・タンパク質中のヘテロ原子座標や側鎖構造の最適化

・MDの前処理(preparemd by Moriwaki Yoshitaka san)

・jobの投入(qsub.sh)や出力結果の取得(init.sh)を一括で行うスクリプトの作成

・実験条件のバックアップ(config.ini)

を、複数のIDを対象にワンストップで行うスクリプト

## Requirements
以下のパッケージを実行可能なPython3.9.x 環境(Pyrosettaが3.9のみでしかインポートできないため)

・metapredict

・Biopython

・Pymol2

・Pyrosetta

・absl

・preparemd

・modeller

## How to Use
以下のコマンドで実行。

    python3.9 ./pdb2md.py --config_file ./config.ini

引数は1つで、--config_fileでconfig.iniのパスを指定する。デフォルト値は./config.ini

## Setting up config.ini

処理は次の順番で行われる
    
    ・指定したIDの結晶構造 もしくはAlphaFoldによる予測構造を取得
    ↓
    ・*(結晶構造に対して)ミッシング領域の構築
    ↓
    ・*ディスオーダー領域の削除
    ↓
    ・*テンプレート構造を基にヘテロ原子を挿入
    ↓
    ・*ヘテロ原子座標や側鎖構造の最適化
    ↓
    ・*テンプレート構造基に基質を挿入
    ↓
    ・MDの前処理(preparemd by Moriwaki Yoshitaka san)
    ↓
    ・jobの投入(qsub.sh)や出力結果の取得(init.sh)を一括で行うスクリプトの作成&実験条件のバックアップ(config.ini)
    
    * がついた処理はconfig.iniからon/offを切り替えられる

基本的にconfig.iniを編集して操作する。凡例はconfig_templete.iniを参照。

### [PATH]

    distination_path = 作業フォルダを作成する場所

    templete_pdb_path = テンプレートとするpdbファイルの場所

    preparemd_script_path = preparemdのスクリプトの場所

    parameter_file_path = 基質のパラメータファイルの場所

### [SETTINGS]

    workbench_dir_name = 作成する作業フォルダの名前

    remove_disordered_residue = ディスオーダー領域を削除するか

    using_modeller_for_disordered_residue = modellerで構造内部のミッシング領域を構築するかどうか

    insert_residue_from_temolete = テンプレートのヘテロ原子などを挿入するかどうか
    ※具体的な原子もしくは残基名は次セクション[RESIDUES_NAME_IN_TEMPLETE]で指定する

    rosetta_packing = Rosettaで基質・ヘテロ金属原子の最適化を行うかどうか

    insert_substrate_from_templete = テンプレートの基質を挿入するかどうか
    ※具体的な基質名は次セクション[RESIDUES_NAME_IN_TEMPLETE]で指定する

### [RESIDUES_NAME_IN_TEMPLETE]

    insert_residue_name = 挿入するヘテロ原子などの名前(複数種未対応)

    insert_substrate_name = 挿入する基質の名前(複数未対応)

### [ROSETTA_SETTINGS]

 Rosettaのinit()で指定する引数 (「-」「--」は除いて指定すること)
 オプションが引数を取る場合は=で指定(ex ignore_zero_occupancy = False)

### [PREPAREMD_SETTINGS]

 preparemdの引数。<https://github.com/YoshitakaMo/preparemd> に準拠しているためそちらを参照
 (「-」「--」は除いて指定すること)

### [ID]

対象とするuniplot/PDB IDを指定する。

複数指定可能。

大文字小文字のどちらでも指定可能。

重複したIDは片方のみ処理される。

## Future Work

## Tips
