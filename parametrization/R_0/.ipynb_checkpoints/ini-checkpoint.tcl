mol delete all
mol new ../template_R.pdb
set all [atomselect top "not index 20"]
$all set resname RN
$all writepdb SM86-R-0.pdb
exit
