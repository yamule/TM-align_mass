# TM-align_mass
TM-align https://zhanggroup.org/TM-align/ for large dbs.

## Main differences with the original TM-align:
 - Check alignment state in one more previous step in TMscore8_search(_standard) function.
 - Force one of the NWDP_TM function to be inline.
 - Added binary file mode.

Anyway, I expect that result is completely the same with that of the original TM-align.

The above modifications reduced processing time just 8% or so, therefore, you can use the original version with this small loss of processing time.

## Compile:
```
clone https://github.com/yamule/TM-align_mass.git
cd TM-align_mass
mkdir bin
g++ -static -O3 -ffast-math -lm -o bin/TMalign_mass.exe TMalign_mass.cpp -lz
```

## An example usage:
```
mkdir UP000000625_83333_ECOLI
cd UP000000625_83333_ECOLI
wget https://ftp.ebi.ac.uk/pub/databases/alphafold/UP000000625_83333_ECOLI.tar
tar xvf UP000000625_83333_ECOLI.tar
gunzip *.pdb.gz
ls *.pdb > pdblist.dat
cd ..
bin/TMalign_mass.exe UP000000625_83333_ECOLI/AF-P36677-F1-model_v1.pdb -dir2 UP000000625_83333_ECOLI/ UP000000625_83333_ECOLI/pdblist.dat
```


## Perform multi-processing:
(You have to be in the directory of this repo & finished "An example usage" except for the last line.)
```
find UP000000625_83333_ECOLI | grep -E "\.pdb$" | xargs -I {} bin/TMalign_mass.exe {} -bin_convert {}.bxyz 
ls UP000000625_83333_ECOLI | grep -E "\.pdb.bxyz$"  > UP000000625_83333_ECOLI/pdblist_binary.dat
perl scripts/multi_process.pl --query UP000000625_83333_ECOLI/AF-P36677-F1-model_v1.pdb --dir UP000000625_83333_ECOLI/ --list UP000000625_83333_ECOLI/pdblist_binary.dat --num_threads 6 --tmscore_threshold 0.5 
```

## Generate clustered db:
(You have to be in the directory of this repo & finished "Perform multi-processing" except for the last line.)

Set path for blastp and makeblastdb in scripts/iter_blastclust.pl
```
mkdir output
perl scripts/getseq_pdb.pl UP000000625_83333_ECOLI/ > output/allfas.fas
perl scripts/iter_blastclust.pl -in output/allfas.fas -out output/clustered.dat -outdir output/clusout -ident 30 -cov_long 0.5 -cov_short 0.5
ここから
クラスタ内の構造比較
一時ファイルはTmpDirに作るようにする
```


## References to cite:
 - TM-align<br>
Y Zhang, J Skolnick. Nucl Acids Res 33, 2302-9 (2005)
https://academic.oup.com/nar/article/33/7/2302/2401364

 - This repo<br>
https://github.com/yamule/TM-align_mass

## License:
 - TM-align<br>
Please follow the comment in the head of the source code.
 - Scripts under the scripts/ dir<br>
Apache License, Version 2.0
http://www.apache.org/licenses/LICENSE-2.0.html
