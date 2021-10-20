# TM-align_mass
TM-align for large dbs.

## Compile:
```
clone https://github.com/yamule/TM-align_mass.git
cd TM-align_mass
mkdir bin
g++ -static -O3 -ffast-math -lm -o bin/TMalign_mass.exe TMalign_mass.cpp
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
bin/TMalign_mass.exe UP000000625_83333_ECOLI/AF-P38097-F1-model_v1.pdb -dir2 UP000000625_83333_ECOLI/ UP000000625_83333_ECOLI/pdblist.dat
```

## References to cite:
Y Zhang, J Skolnick. Nucl Acids Res 33, 2302-9 (2005)


## License:
Please follow the comment in the head of the source code.
