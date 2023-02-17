#Centering lipids in the simulation box
#The system has to be built by charmgui and for NAMD simulation
proc centering {} {
	set sel [atomselect top "segname MEMB"] 
	set all [atomselect top "all"]
	pbc wrap -center bb -centersel $sel -compound res -all
	set r1 [measure center $sel]
	set r2 [vecinvert $r1]
	$all moveby $r2
}
#resids of lipids in the upper leaflet chosen to be exchanged
#if the novel lipid is not phospholipid, "name P" needs to change to something else
proc lipid_resids_upper {exchange_selectin number_novel_upper} {
	set upper_exchange_selectin [atomselect top "resname $exchange_selectin and name P and z > 0"]
	set resids_upper [$upper_exchange_selectin get resid]
	set N1 [$upper_exchange_selectin num]
	set List_upper {}

	for {set i 0} {$i < $number_novel_upper} {incr i} {
		set j [expr {int(rand() * $N1)}]
		lappend List_upper [lindex $resids_upper $j]
	}
	return $List_upper
}
#resids of lipids in the lower leaflet chosen to be exchanged
proc lipid_resids_lower {exchange_selectin number_novel_lower} {

	set lower_exchange_selectin [atomselect top "resname $exchange_selectin and name P and z < 0"]
	set resids_lower [$lower_exchange_selectin get resid]
	set N1 [$lower_exchange_selectin num]
	set List_lower {}

	for {set i 0} {$i < $number_novel_lower} {incr i} {
		set j [expr {int(rand() * $N1)}]
		lappend List_lower [lindex $resids_lower $j]
	}
	return $List_lower
}
#COM of each lipid in the upper leaflet chosen to be exchanged
proc lipids_COM_upper {exchange_selectin lipid_resids_upper} {
	set COM_upper {}
	foreach resid $lipid_resids_upper {
		set sel [atomselect top "resname $exchange_selectin and resid $resid"]
		set COM [measure center $sel]
		lappend COM_upper $COM
	}
	return $COM_upper
}
#COM of each lipid in the lower leaflet chosen to be exchanged
proc lipids_COM_lower {exchange_selectin lipid_resids_lower} {
	set COM_lower {}
	foreach resid $lipid_resids_lower {
		set sel [atomselect top "resname $exchange_selectin and resid $resid"]
		set COM [measure center $sel]
		lappend COM_lower $COM
	}
	return $COM_lower
}
#Drwing a principal axes for a selection
proc drawprincipalaxes {R a b} {
    set CM [lindex $R 0]
    set R1 [vecscale $a [lindex $R 1 0]]
    set R2 [vecscale $a [lindex $R 1 1]]
    set R3 [vecscale $a [lindex $R 1 2]]
    set R4 [vecscale $b [lindex $R 1 0]]
    set R5 [vecscale $b [lindex $R 1 1]]
    set R6 [vecscale $b [lindex $R 1 2]]

    set CM1 [vecadd $CM $R1]
    set CM2 [vecadd $CM $R2]
    set CM3 [vecadd $CM $R3]
    set CM4 [vecadd $CM1 $R4]
    set CM5 [vecadd $CM2 $R5]
    set CM6 [vecadd $CM3 $R6]

    draw delete all
    draw color blue2
    draw cylinder $CM $CM1 radius 0.1
    draw cone $CM1 $CM4 radius 0.2
    draw color green2
    draw cylinder $CM $CM2 radius 0.1
    draw cone $CM2 $CM5 radius 0.2
    draw color red
    draw cylinder $CM $CM3 radius 0.1
    draw cone $CM3 $CM6 radius 0.2

}
#preperating templates for lower and upper leaflets
proc template {pdbfile psffile counter} {
	mol load psf $psffile pdb $pdbfile
	set sel [atomselect top all]
	set r1 [measure center $sel]
	set r2 [vecinvert $r1]
	$sel moveby $r2
	#set I [measure inertia $sel]
	# drawprincipalaxes $I 2 2
	#set sel1 [atomselect top "name N"]
	#set sel2 [atomselect top "name O O4"]
	#set COM1 [measure center $sel1]
	#set COM2 [measure center $sel2]
	#set r_ini [vecsub $COM1 $COM2]
	#set M  [transvecinv $r_ini]
	#set M1 [trans axis y 90]
	set M2 [trans axis y +180]

	#$sel move $M
	#$sel move $M1
	$sel writepdb ../templates/template-$counter-upper.pdb

	$sel move $M2
	$sel writepdb ../templates/template-$counter-lower.pdb
	mol delete all
}
#merging structures
proc merging_structure {input output N} {
	mol delete all
	package require psfgen 
	resetpsf

	for {set i 0} {$i < $N} {incr i} {
		readpsf  $input-$i.psf
		coordpdb $input-$i.pdb
	}
	writepsf $output.psf
	writepdb $output.pdb
}
#generating colvar files for phospholipids and cholestrols
#generating colvar files for phospholipids and cholestrols
proc colvars_0 {} {
	set all [atomselect top all]
	$all set beta 0
	$all set occupancy 0
	set sel2_l [atomselect top "segname MEMB and name P and z < 0"]
	set sel2_u [atomselect top "segname MEMB and name P and z > 0"]
	set sel3_l [atomselect top "segname MEMB and resname CHL1 and name O3 and z < 0"]
	set sel3_u [atomselect top "segname MEMB and resname CHL1 and name O3 and z > 0"]

	$all set beta 0
	$all set occupancy 0
	$sel2_l set beta 1
	$all writepdb PH_head_lower_1.pdb

	$all set beta 0
	$all set occupancy 0
	$sel3_l set beta 1
	$all writepdb CH_head_lower_1.pdb

	$all set beta 0
	$all set occupancy 0
	$sel2_u set beta 1
	$all writepdb PH_head_upper_1.pdb

	$all set beta 0
	$all set occupancy 0
	$sel3_u set beta 1
	$all writepdb CH_head_upper_1.pdb

	$all set beta 0
	$all set occupancy 0

	set fp [open "colvar_0.conf" w+]
	puts $fp "###Lipid Head###"
	puts $fp "Colvarstrajfrequency    100"
	puts $fp "Colvarsrestartfrequency 100"
	puts $fp "colvar {"
   	puts $fp "name chl1_head_upper"
   	puts $fp "distanceZ {"
   	puts $fp "ref {"
   	puts $fp "dummyAtom ( 0.000, 0.000, 0.000 )"
   	puts $fp "}"
   	puts $fp "main {"
   	puts $fp "atomsFile      CH_head_upper_1.pdb"
   	puts $fp "atomsCol       B"
   	puts $fp "atomsColValue  1.0"
   	puts $fp "}"
   	puts $fp "}"
   	puts $fp "}"
	puts $fp "harmonic {"
	puts $fp "colvars chl1_head_upper"
	puts $fp "centers 18"
	puts $fp "forceConstant 5"
	puts $fp "}"
	puts $fp "colvar {"
	puts $fp "name phos_head_upper"
	puts $fp "distanceZ {"
	puts $fp "ref {"
	puts $fp "dummyAtom ( 0.000, 0.000, 0.000 )"
	puts $fp "}"
	puts $fp "main {"
	puts $fp "atomsFile      PH_head_upper_1.pdb"
	puts $fp "atomsCol       B"
	puts $fp "atomsColValue  1.0"
	puts $fp "}"
	puts $fp "}"
	puts $fp "}"
	puts $fp "harmonic {"
	puts $fp "colvars phos_head_upper"
	puts $fp "centers 19"
	puts $fp "forceConstant 5"
	puts $fp "}"
	puts $fp "colvar {"
   	puts $fp "name chl1_head_lower"
   	puts $fp "distanceZ {"
   	puts $fp "ref {"
   	puts $fp "dummyAtom ( 0.000, 0.000, 0.000 )"
   	puts $fp "}"
   	puts $fp "main {"
   	puts $fp "atomsFile      CH_head_lower_1.pdb"
   	puts $fp "atomsCol       B"
   	puts $fp "atomsColValue  1.0"
   	puts $fp "}"
   	puts $fp "}"
   	puts $fp "}"
	puts $fp "harmonic {"
	puts $fp "colvars chl1_head_lower"
	puts $fp "centers -18"
	puts $fp "forceConstant 5"
	puts $fp "}"
	puts $fp "colvar {"
	puts $fp "name phos_head_lower"
	puts $fp "distanceZ {"
	puts $fp "ref {"
	puts $fp "dummyAtom ( 0.000, 0.000, 0.000 )"
	puts $fp "}"
	puts $fp "main {"
	puts $fp "atomsFile      PH_head_lower_1.pdb"
	puts $fp "atomsCol       B"
	puts $fp "atomsColValue  1.0"
	puts $fp "}"
	puts $fp "}"
	puts $fp "}"
	puts $fp "harmonic {"
	puts $fp "colvars phos_head_lower"
	puts $fp "centers -19"
	puts $fp "forceConstant 5"
	puts $fp "}"
	close $fp          
}


proc conf_0 {toppar pdb psf restart FC} {
	
	set fp [open "gsmd.0.conf" w+]


	puts $fp "structure $psf"
	puts $fp "coordinates $pdb"     
	puts $fp "firsttimestep   0"


	puts $fp "binCoordinates $restart.restart.coor"
	puts $fp "binVelocities  $restart.restart.vel"
	puts $fp "extendedSystem $restart.restart.xsc"


	puts $fp "dcdfreq        100"
	puts $fp "outputEnergies 100"
	puts $fp "xstFreq        100"
	puts $fp "restartfreq    100"
	puts $fp "outputTiming   100"
	puts $fp "outputName     gsmd.0"

	puts $fp "paraTypeCharmm     on"

	puts $fp "parameters     $toppar/par_all36_lipid.prm"
	puts $fp "parameters     $toppar/par_all36_cgenff.prm"
	puts $fp "parameters     $toppar/toppar_water_ions.str"
	puts $fp "parameters     $toppar/toppar_all36_lipid_bacterial.str"
	puts $fp "parameters     $toppar/toppar_all36_lipid_cardiolipin.str"
	puts $fp "parameters     $toppar/toppar_all36_lipid_cholesterol.str"
	puts $fp "parameters     $toppar/toppar_all36_lipid_inositol.str"
	puts $fp "parameters     $toppar/toppar_all36_lipid_lps.str"
	puts $fp "parameters     $toppar/toppar_all36_lipid_miscellaneous.str"
	puts $fp "parameters     $toppar/toppar_all36_lipid_model.str"
	puts $fp "parameters     $toppar/toppar_all36_lipid_prot.str"
	puts $fp "parameters     $toppar/toppar_all36_lipid_sphingo.str"
	puts $fp "parameters     $toppar/toppar_all36_lipid_yeast.str"
	puts $fp "parameters     $toppar/toppar_all36_lipid_hmmm.str"
	puts $fp "parameters     $toppar/toppar_all36_lipid_detergent.str"
	puts $fp "parameters     $toppar/toppar_all36_lipid_ether.str"
	puts $fp "parameters     $toppar/par_all36m_prot.prm"
	puts $fp "parameters     $toppar/par_all36_na.prm"
	puts $fp "parameters     $toppar/par_all36_carb.prm"
	puts $fp "parameters     $toppar/par_interface.prm"
	puts $fp "parameters     $toppar/toppar_all36_nano_lig.str"
	puts $fp "parameters     $toppar/toppar_all36_nanolig_patch.str"
	puts $fp "parameters     $toppar/toppar_dum_noble_gases.str"
	puts $fp "parameters     $toppar/toppar_ions_won.str"
	puts $fp "parameters     $toppar/toppar_all36_prot_c36m_d_aminoacids.str"
	puts $fp "parameters     $toppar/toppar_all36_prot_fluoro_alkanes.str"
	puts $fp "parameters     $toppar/toppar_all36_prot_na_combined.str"
	puts $fp "parameters     $toppar/toppar_all36_prot_heme.str"
	puts $fp "parameters     $toppar/toppar_all36_prot_retinol.str"
	puts $fp "parameters     $toppar/toppar_all36_prot_modify_res.str"
	puts $fp "parameters     $toppar/toppar_all36_na_nad_ppi.str"
	puts $fp "parameters     $toppar/toppar_all36_na_rna_modified.str"
	puts $fp "parameters     $toppar/toppar_all36_carb_glycolipid.str"
	puts $fp "parameters     $toppar/toppar_all36_carb_glycopeptide.str"
	puts $fp "parameters     $toppar/toppar_all36_carb_imlab.str"
	puts $fp "parameters     $toppar/toppar_all36_label_spin.str"
	puts $fp "parameters     $toppar/toppar_all36_label_fluorophore.str"

	puts $fp "exclude             scaled1-4"
	puts $fp "1-4scaling          1.0"
	puts $fp "switching           on"
	puts $fp "vdwForceSwitching   yes"

	puts $fp "cutoff              12.0"
	puts $fp "switchdist          10.0"
	puts $fp "pairlistdist        13.5"
	puts $fp "stepspercycle       20"
	puts $fp "pairlistsPerCycle    2"
                          
	puts $fp "timestep            2.0"
	puts $fp "rigidBonds          all"
	puts $fp "nonbondedFreq       1"
	puts $fp "fullElectFrequency  2"

	puts $fp "wrapWater           on"
	puts $fp "wrapAll             on"

	puts $fp "PME                yes"
	puts $fp "PMEInterpOrder       8"
	puts $fp "PMEGridSpacing     2.1"
	puts $fp "PMEOffload 	   yes"

	puts $fp "useGroupPressure       yes"
	puts $fp "useFlexibleCell        yes"
	puts $fp "useConstantRatio       yes"

	puts $fp "langevinPiston          on"
	puts $fp "langevinPistonTarget  1.01325"
	puts $fp "langevinPistonPeriod  200.0"
	puts $fp "langevinPistonDecay   200.0"
	puts $fp "langevinPistonTemp   310"

	puts $fp "langevin                on"
	puts $fp "langevinDamping        1.0"
	puts $fp "langevinTemp         310"
	puts $fp "langevinHydrogen       off"

	puts $fp "gridforce              yes"
	puts $fp "gridforcefile          ref_0.pdb"
	puts $fp "gridforcecol           B"
	puts $fp "gridforcechargecol     O"
	puts $fp "gridforcepotfile       volmap_out_0.dx"
	puts $fp "gridforcescale         $FC $FC $FC"
	puts $fp  "gridforcechecksize    off"

	puts $fp "colvars                 on"
	puts $fp "colvarsConfig           colvar_0.conf"

	puts $fp "run 1000"  
	close $fp          
}


proc conf_2 {filename pdb par restart FC1 outputname} {
	
	set fp [open "$filename" w+]


	puts $fp "parmfile $par"
	puts $fp "coordinates $pdb"     
	puts $fp "firsttimestep   0"

        puts $fp "binCoordinates $restart.restart.coor"
        puts $fp "binVelocities  $restart.restart.vel"
        puts $fp "extendedSystem $restart.restart.xsc"
	
	puts $fp "dcdfreq        100"
	puts $fp "outputEnergies 100"
	puts $fp "xstFreq        100"
	puts $fp "restartfreq    100"
	puts $fp "outputTiming   100"
	puts $fp "outputName     $outputname"

	puts $fp "amber     on"
	puts $fp "wrapWater on"
	puts $fp "wrapAll on"
	
	puts $fp "exclude             scaled1-4"
	puts $fp "1-4scaling          0.833333"
	puts $fp "cutoff              9.0"
	puts $fp "switching off"
	puts $fp "pairlistdist        11.0"
	puts $fp "readexclusions      yes"
	puts $fp "stepspercycle       10"
	puts $fp "scnb 2.0"
	puts $fp "watermodel tip3"
	puts $fp "LJcorrection on"
	puts $fp "ZeroMomentum on"
                   
	puts $fp "timestep            2.0"
	puts $fp "rigidBonds          all"
	puts $fp "rigidTolerance  1.0e-8"
	puts $fp "rigidIterations 100"
	puts $fp "UseSettle on"

	puts $fp "nonbondedFreq       1"
	puts $fp "fullElectFrequency  1"

	puts $fp "PME                yes"
	puts $fp "PMEInterpOrder       6"
	puts $fp "PMEGridSpacing     1.0"
	puts $fp "PMETolerance       1.0e-6"

	puts $fp "useGroupPressure       yes"
	puts $fp "useFlexibleCell        yes"
	puts $fp "useConstantRatio       yes"

	puts $fp "langevinPiston          on"
	puts $fp "langevinPistonTarget  1.01325"
	puts $fp "langevinPistonPeriod  100.0"
	puts $fp "langevinPistonDecay   50.0"
	puts $fp "langevinPistonTemp   310"

	puts $fp "langevin                on"
	puts $fp "langevinDamping        1.0"
	puts $fp "langevinTemp         310"
	puts $fp "langevinHydrogen       off"

	puts $fp "gridforce              yes"
	puts $fp "gridforcefile          ref_1.pdb"
	puts $fp "gridforcecol           B"
	puts $fp "gridforcechargecol     O"
	puts $fp "gridforcepotfile       volmap_out_0.dx"
	puts $fp "gridforcescale         5.0 5.0 5.0"
	puts $fp  "gridforcechecksize    off"

	puts $fp "colvars                 on"
	puts $fp "colvarsConfig           colvar_1.conf"

	puts $fp "updateGridScale $FC1 $FC1 $FC1"
	puts $fp "run 1000"  

	close $fp          
}


proc conf_1 {pdb par restart FC1 r1 r2 r3} {
	
	set fp [open "gsmd.1.conf" w+]


	puts $fp "parmfile $par"
	puts $fp "coordinates $pdb"     
	puts $fp "firsttimestep   0"

	puts $fp "cellBasisVector1     $r1   0.0   0.0"
	puts $fp "cellBasisVector2     0.0   $r2   0.0"
	puts $fp "cellBasisVector3     0.0   0.0   $r3"
	puts $fp "cellOrigin           0.0   0.0   0.0"
#	puts $fp "extendedSystem $restart.restart.xsc"
	puts $fp "margin 10"
        puts $fp "temperature 310"

	puts $fp "dcdfreq        100"
	puts $fp "outputEnergies 100"
	puts $fp "xstFreq        100"
	puts $fp "restartfreq    100"
	puts $fp "outputTiming   100"
	puts $fp "outputName     gsmd.1"

	puts $fp "amber     on"
	puts $fp "wrapWater on"
	puts $fp "wrapAll on"
	
	puts $fp "exclude             scaled1-4"
	puts $fp "1-4scaling          0.833333"
	puts $fp "cutoff              9.0"
	puts $fp "switching off"
	puts $fp "pairlistdist        11.0"
	puts $fp "readexclusions      yes"
	puts $fp "stepspercycle       10"
	puts $fp "scnb 2.0"
	puts $fp "watermodel tip3"
	puts $fp "LJcorrection on"
	puts $fp "ZeroMomentum on"
                   
	puts $fp "timestep            2.0"
	puts $fp "rigidBonds          all"
	puts $fp "rigidTolerance  1.0e-8"
	puts $fp "rigidIterations 100"
	puts $fp "UseSettle on"

	puts $fp "nonbondedFreq       1"
	puts $fp "fullElectFrequency  1"

	puts $fp "PME                yes"
	puts $fp "PMEInterpOrder       6"
	puts $fp "PMEGridSpacing     1.0"
	puts $fp "PMETolerance       1.0e-6"

	puts $fp "useGroupPressure       yes"
	puts $fp "useFlexibleCell        yes"
	puts $fp "useConstantRatio       yes"

	puts $fp "langevinPiston          on"
	puts $fp "langevinPistonTarget  1.01325"
	puts $fp "langevinPistonPeriod  100.0"
	puts $fp "langevinPistonDecay   50.0"
	puts $fp "langevinPistonTemp   310"

	puts $fp "langevin                on"
	puts $fp "langevinDamping        1.0"
	puts $fp "langevinTemp         310"
	puts $fp "langevinHydrogen       off"

	puts $fp "gridforce              yes"
	puts $fp "gridforcefile          ref_1.pdb"
	puts $fp "gridforcecol           B"
	puts $fp "gridforcechargecol     O"
	puts $fp "gridforcepotfile       volmap_out_0.dx"
	puts $fp "gridforcescale         5.0 5.0 5.0"
	puts $fp  "gridforcechecksize    off"

	puts $fp "colvars                 on"
	puts $fp "colvarsConfig           colvar_1.conf"

	puts $fp "minimize 1000"
#	puts $fp "updateGridScale $FC1 $FC1 $FC1"
#	puts $fp "run 1000"  

	close $fp          
}


proc colvars_1 {} {
	set all [atomselect top all]
	$all set beta 0
	$all set occupancy 0
        
	set sel2_l [atomselect top "resname PC and name P31 and z < 0"]
        set sel2_u [atomselect top "resname PC and name P31 and z > 0"]
        set sel3_l [atomselect top "resname CHL and name O1 and z < 0"]
        set sel3_u [atomselect top "resname CHL and name O1 and z > 0"]

	$all set beta 0
	$all set occupancy 0
	$sel2_l set beta 1
	$all writepdb PH_head_lower_2.pdb

	$all set beta 0
	$all set occupancy 0
	$sel3_l set beta 1
	$all writepdb CH_head_lower_2.pdb

	$all set beta 0
	$all set occupancy 0
	$sel2_u set beta 1
	$all writepdb PH_head_upper_2.pdb

	$all set beta 0
	$all set occupancy 0
	$sel3_u set beta 1
	$all writepdb CH_head_upper_2.pdb

	$all set beta 0
	$all set occupancy 0

	set fp [open "colvar_1.conf" w+]
	puts $fp "###Lipid Head###"
	puts $fp "Colvarstrajfrequency    100"
	puts $fp "Colvarsrestartfrequency 100"
	puts $fp "colvar {"
   	puts $fp "name chl1_head_upper"
   	puts $fp "distanceZ {"
   	puts $fp "ref {"
   	puts $fp "dummyAtom ( 0.000, 0.000, 0.000 )"
   	puts $fp "}"
   	puts $fp "main {"
   	puts $fp "atomsFile      CH_head_upper_2.pdb"
   	puts $fp "atomsCol       B"
   	puts $fp "atomsColValue  1.0"
   	puts $fp "}"
   	puts $fp "}"
   	puts $fp "}"
	puts $fp "harmonic {"
	puts $fp "colvars chl1_head_upper"
	puts $fp "centers 18"
	puts $fp "forceConstant 5"
	puts $fp "}"
	puts $fp "colvar {"
	puts $fp "name phos_head_upper"
	puts $fp "distanceZ {"
	puts $fp "ref {"
	puts $fp "dummyAtom ( 0.000, 0.000, 0.000 )"
	puts $fp "}"
	puts $fp "main {"
	puts $fp "atomsFile      PH_head_upper_2.pdb"
	puts $fp "atomsCol       B"
	puts $fp "atomsColValue  1.0"
	puts $fp "}"
	puts $fp "}"
	puts $fp "}"
	puts $fp "harmonic {"
	puts $fp "colvars phos_head_upper"
	puts $fp "centers 19"
	puts $fp "forceConstant 5"
	puts $fp "}"
	puts $fp "colvar {"
   	puts $fp "name chl1_head_lower"
   	puts $fp "distanceZ {"
   	puts $fp "ref {"
   	puts $fp "dummyAtom ( 0.000, 0.000, 0.000 )"
   	puts $fp "}"
   	puts $fp "main {"
   	puts $fp "atomsFile      CH_head_lower_2.pdb"
   	puts $fp "atomsCol       B"
   	puts $fp "atomsColValue  1.0"
   	puts $fp "}"
   	puts $fp "}"
   	puts $fp "}"
	puts $fp "harmonic {"
	puts $fp "colvars chl1_head_lower"
	puts $fp "centers -18"
	puts $fp "forceConstant 5"
	puts $fp "}"
	puts $fp "colvar {"
	puts $fp "name phos_head_lower"
	puts $fp "distanceZ {"
	puts $fp "ref {"
	puts $fp "dummyAtom ( 0.000, 0.000, 0.000 )"
	puts $fp "}"
	puts $fp "main {"
	puts $fp "atomsFile      PH_head_lower_2.pdb"
	puts $fp "atomsCol       B"
	puts $fp "atomsColValue  1.0"
	puts $fp "}"
	puts $fp "}"
	puts $fp "}"
	puts $fp "harmonic {"
	puts $fp "colvars phos_head_lower"
	puts $fp "centers -19"
	puts $fp "forceConstant 5"
	puts $fp "}"
	close $fp          
}

proc generate_psf {pdbfile parfile} {

        set fp [open "par2psf.py" w+]
	puts $fp "import parmed as pmd"
	puts $fp "amber = pmd.load_file('$parfile', '$pdbfile')"
	puts $fp "amber.save('initial_system.psf')"

	close $fp
}

proc generate_par {pdbfile psffile} {

        set fp [open "par2psf.py" w+]
        puts $fp "import parmed as pmd"
        puts $fp "amber = pmd.load_file('$psffile', '$pdbfile')"
        puts $fp "amber.save('initial_system.parm7')"

        close $fp
}

proc tleap {} {
	set fp [open "tleap.in" w+]
	puts $fp "source oldff/leaprc.ff99SB"
	puts $fp "source leaprc.gaff"
	puts $fp "source leaprc.lipid17"
	puts $fp "loadoff ../../../parametrization/R_1/RP.lib"
	puts $fp "loadoff ../../../parametrization/S_1/SP.lib"
	puts $fp "loadoff ../../../parametrization/R_0/RN.lib"
	puts $fp "loadoff ../../../parametrization/S_0/SN.lib"
	puts $fp "lipids = loadpdb ionized_vmd_fixed.pdb"
	#puts $fp "solvatebox lipids TIP3PBOX {0 0 15} 1.5"
	#puts $fp "set lipids box {$r1 $r2 $r3}"
	#puts $fp " addionsRand lipids Cl- 0"
	puts $fp "saveamberparm lipids ionized.prmtop ionized.inpcrd"
	puts $fp "quit"

	close $fp
}

proc run_tleap {} {
	set fp [open "run_tleap.sh" w+]
	puts $fp "#!/bin/bash"
	puts $fp "tleap -f tleap.in"

	close $fp
}

proc solvate_ionized_vmd {r1 r2 r3} {
	package require solvate
	set m_z [expr {$r3 / 2.0} - 20]
	solvate new_system.psf new_system.pdb -o solvated -s WF -x 0 -y 0 -z $m_z +x 0 +y 0 +z $m_z -b 2

	file delete combine.pdb
	file delete combine.psf
	file delete solvated.log

	mol delete all

	mol load psf solvated.psf pdb solvated.pdb
	set LIPID [atomselect top "segname MEMB"]
	set v [measure minmax $LIPID ]
	set MAX_Z_L [lindex $v {1 2}]
	set MIN_Z_L [lindex $v {0 2}]
#	set MAX_X [lindex $v {1 0}]
#	set MIN_X [lindex $v {0 0}]
#	set MAX_Y [lindex $v {1 1}]
#	set MIN_Y [lindex $v {0 1}]
	set MAX_Z [expr {$r3 /  2}]
	set MIN_Z [expr {$r3 / -2}]
        set MAX_Y [expr {$r2 /  2}]
        set MIN_Y [expr {$r2 / -2}]
	set MAX_X [expr {$r1 /  2}]
        set MIN_X [expr {$r1 / -2}]
#	set seltext1 "(same residue as water and (z < $MAX_Z and z > $MIN_Z and x < $MAX_X and x > $MIN_X and y < $MAX_Y and y > $MIN_Y))"
	set seltext1 "(same residue as water and (z < $MAX_Z_L and z > $MIN_Z_L))"
	set seltext2 "(same residue as water and (x > $MAX_X or x < $MIN_X or y > $MAX_Y or y < $MIN_Y))"
	set watSel [atomselect top "not ($seltext1) and not ($seltext2)"]


	$watSel writepdb solvated2.pdb
	$watSel writepsf solvated2.psf

	package require autoionize
        autoionize -psf solvated2.psf -pdb solvated2.pdb -o ionized_vmd_neutral -neutralize
        autoionize -psf ionized_vmd_neutral.psf -pdb ionized_vmd_neutral.pdb -o ionized_vmd -sc 0.150
	::ExecTool::exec solvated2.psf solvated2.pdb solvated.psf solvated.pdb
}
### 


### Removing Ring piercing
proc ring_piercing {cutoff psffile dcdfile coordfile} {

	set fp [open "$psffile" r]
	set lines [split [read $fp] "\n"]
	close $fp
	set N [llength $lines]

	for {set i 0} {$i < $N} {incr i} {
		set line [lindex $lines $i]
		if {[lindex $line 1] == "!NBOND:"} {
			set I1 $i
		}
		if {[lindex $line 1] == "!NTHETA:"} {
			set I2  $i
		}
	}
	set I1 [expr {$I1 + 1}]
	set I2 [expr {$I2 - 1}]

	set list1 {}
	set list2 {}
	
	mol delete all
	mol load psf $psffile
	mol addfile $dcdfile waitfor all
	set all [atomselect top all]
	for {set i $I1} {$i < $I2} {incr i} {

		set line [lindex $lines $i]
		set L    [expr {[llength $line] * 0.5}]

		for {set j 0} {$j < $L} {incr j} {

			set J1 [expr {$j * 2}]
			set J2 [expr {$j * 2 + 1}]

			set first   [expr 1 * [lindex $line $J1] - 1]
			set second  [expr 1 * [lindex $line $J2] - 1]
			# puts $first 
			# puts $second

			set r [measure bond "$first $second"]
			if {$r > $cutoff} {
				lappend list1 $first
				lappend list2 $second
				puts "there is a ring piercing :-(, atoms $first and $second are involved"
			}
		}
	}
	foreach first $list1 second $list2 {
		set sel1 [atomselect top "index $first"]
		set sel2 [atomselect top "index $second"]
		set COM1 [measure center $sel1]
		set COM2 [measure center $sel2]
		set vec1 [vecsub $COM1 $COM2]
		set x 	 [vecinvert [lindex $vec1 1]]
		set y 	 [lindex $vec1 0]
		set vect [list $x $y 0]
		$sel1 moveby $vect
		$sel2 moveby $vect
	}
	$all writenamdbin $coordfile

}

proc generate_psf {pdbfile parfile} {

        set fp [open "par2psf.py" w+]
        puts $fp "import parmed as pmd"
        puts $fp "amber = pmd.load_file('$parfile', '$pdbfile')"
        puts $fp "amber.save('ionized.psf')"

        close $fp
}
