 #!/bin/bash
#conda activate AmberTools20
rm ANTECHAMBER* *pdb *psf *lib sqm* *log *INF *prmtop *inpcrd *mol2 *frcmod
vmd -dispdev text -e ini.tcl
antechamber -i SM86-R-0.pdb  -fi pdb -o SM86-R-0.mol2   -fo mol2 -c bcc -nc 0 -s 2
parmchk2    -i SM86-R-0.mol2 -f mol2 -o SM86-R-0.frcmod
tleap -f tleap.in
