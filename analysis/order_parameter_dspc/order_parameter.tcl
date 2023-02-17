
mol delete all
set Files1 [list "10"  "20"    "30"  "40"    "50"   ]
set Files2 [list "1"   "1"     "1"   "1"     "1"    ]
set Files3 [list "equ" "equ_1" "equ" "equ_1" "equ_1"]
foreach f1 $Files1 f2 $Files2 f3 $Files3 {
    mol load psf "../../../$f1/$f2/output/ionized.psf"
    mol addfile "../../../$f1/$f2/$f3/reduced.dcd" waitfor all step 1
    set output0 [open order_parameter_$f1.txt w]

    set nf [molinfo top get numframes]
    set dspc_selection [atomselect top "resname PC and name P31"]
    set dspc_resid [$dspc_selection get resid]    

    set start [expr {int($nf * 0.666666666666666666)}]
    
    for {set f $start} {$f <$nf} {incr f} {
    
        foreach resid $dspc_resid {
            set resid1 [expr {$resid + 1}]
            set resid2 [expr {$resid - 1}]
            
            set sel1 [atomselect top "resname PC and resid $resid"  frame $f]
            set sel2 [atomselect top "resname ST and resid $resid1 $resid2" frame $f]
            
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


        puts "$f out of $nf"
    }

    close $output0
}
