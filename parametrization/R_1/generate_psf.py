import parmed as pmd

amber = pmd.load_file('SM86-R-1.prmtop', 'SM86-R-1.pdb')

# Save a CHARMM PSF and crd file
amber.save('SM86-R-1.psf')
