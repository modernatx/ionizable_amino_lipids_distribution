 #!/bin/bash
#conda activate AmberTools20
rm ANTECHAMBER* *pdb *psf *lib sqm* *log *INF *prmtop *inpcrd *mol2 *frcmod
vmd -dispdev text -e ini.tcl
antechamber -i SM86-S-1.pdb  -fi pdb -o SM86-S-1.mol2   -fo mol2 -c bcc -nc +1 -s 2
parmchk2    -i SM86-S-1.mol2 -f mol2 -o SM86-S-1.frcmod
tleap -f tleap.in
