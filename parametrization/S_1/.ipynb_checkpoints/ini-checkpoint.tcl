mol delete all
mol new ../template_S.pdb
set all [atomselect top all]
$all set resname SP
$all writepdb SM86-S-1.pdb
exit
