import parmed as pmd

amber = pmd.load_file('SM86-R-0.prmtop', 'SM86-R-0.pdb')

# Save a CHARMM PSF and crd file
amber.save('SM86-R-0.psf')
