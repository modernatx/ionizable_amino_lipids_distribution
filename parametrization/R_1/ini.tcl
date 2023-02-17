mol delete all
mol new ../template_R.pdb
set all [atomselect top all]
$all set resname RP
$all writepdb SM86-R-1.pdb
exit
