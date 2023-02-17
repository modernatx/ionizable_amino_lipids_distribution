
mol delete all
set Files1 [list "10"  "20"    "30"  "40"    "50"   ]
set Files2 [list "1"   "1"     "1"   "1"     "1"    ]
set Files3 [list "equ" "equ_1" "equ" "equ_1" "equ_1"]

foreach f1 $Files1 f2 $Files2 f3 $Files3 {
    mol load psf "../../../$f1/$f2/output/ionized.psf"
    mol addfile "../../../$f1/$f2/$f3/reduced.dcd" waitfor all step 1
    set output0 [open chol_par_$f1.txt w]

    set nf [molinfo top get numframes]
    set initial [atomselect top "resname CHL and name O1"]
    set N [$initial num]
    
    set start [expr {int($nf * 0.6666666666666666666)}]
    
    for {set f $start} {$f <$nf} {incr f} {
    
        set sel1 [atomselect top "resname CHL and name O1 and within 10 of water" frame $f]
        set n [$sel1 num]
        set CHOL_par [expr {1.0 * $n / $N}]
        puts $output0 "$f $CHOL_par"
        $sel1 delete
        
        puts "$f out of $nf"
    }

    close $output0

}
