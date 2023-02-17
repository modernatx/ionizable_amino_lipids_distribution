
mol delete all
set Files1 [list "10"  "20"    "30"  "40"    "50"   ]
set Files2 [list "1"   "1"     "1"   "2"     "1"    ]
set Files3 [list "equ" "equ_1" "equ" "equ_1" "equ_1"]
foreach f1 $Files1 f2 $Files2 f3 $Files3 {
    mol load psf "../../$f1/$f2/output/ionized.psf"
    mol addfile "../../$f1/$f2/$f3/reduced.dcd" waitfor all step 1
    set output0 [open order_parameter_$f1.txt w]
    set output1 [open order_parameter_$f1-p.txt w]
    set output2 [open order_parameter_$f1-n.txt w]

    set nf [molinfo top get numframes]
    set ionizable_selection [atomselect top "resname RN SN RP SP and name N"]
    set ionizable_resid [$ionizable_selection get resid]
    set ionizable_resname [$ionizable_selection get resname]
    

    set ionizable_selection_p [atomselect top "resname RP SP and name N"]
    set ionizable_resid_p [$ionizable_selection_p get resid]
    set ionizable_resname_p [$ionizable_selection_p get resname]


    set ionizable_selection_n [atomselect top "resname RN SN and name N"]
    set ionizable_resid_n [$ionizable_selection_n get resid]
    set ionizable_resname_n [$ionizable_selection_n get resname]

    set start [expr {int($nf * 0.666666666666)}]
    
    for {set f $start} {$f <$nf} {incr f} {
    
        foreach resid $ionizable_resid resname $ionizable_resname {
            set sel1 [atomselect top "resname $resname and resid $resid and name O O4" frame $f]
            set sel2 [atomselect top "resname $resname and resid $resid and name N" frame $f]
            set r1 [measure center $sel1]
            set r2 [measure center $sel2]
            set r3 [vecsub $r2 $r1]
            set rz [lindex $r3 2]
            set rl [veclength $r3]
            set order_parameter [expr {(($rz / $rl) * ($rz / $rl) * 3.0) - 1}]
            puts $output0 "$order_parameter"
            $sel1 delete
            $sel2 delete
        }
        
        foreach resid $ionizable_resid_p resname $ionizable_resname_p {
            set sel1 [atomselect top "resname $resname and resid $resid and name O O4" frame $f]
            set sel2 [atomselect top "resname $resname and resid $resid and name N" frame $f]
            set r1 [measure center $sel1]
            set r2 [measure center $sel2]
            set r3 [vecsub $r2 $r1]
            set rz [lindex $r3 2]
            set rl [veclength $r3]
            set order_parameter [expr {(($rz / $rl) * ($rz / $rl) * 3.0) - 1}]
            puts $output1 "$order_parameter"
            $sel1 delete
            $sel2 delete
        }


        foreach resid $ionizable_resid_n resname $ionizable_resname_n {
            set sel1 [atomselect top "resname $resname and resid $resid and name O O4" frame $f]
            set sel2 [atomselect top "resname $resname and resid $resid and name N" frame $f]
            set r1 [measure center $sel1]
            set r2 [measure center $sel2]
            set r3 [vecsub $r2 $r1]
            set rz [lindex $r3 2]
            set rl [veclength $r3]
            set order_parameter [expr {(($rz / $rl) * ($rz / $rl) * 3.0) - 1}]
            puts $output2 "$order_parameter"
            $sel1 delete
            $sel2 delete
        }

        puts "$f out of $nf"
    }

    close $output0
    close $output1
    close $output2

}
