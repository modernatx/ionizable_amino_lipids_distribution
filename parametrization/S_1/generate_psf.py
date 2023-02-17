import parmed as pmd

amber = pmd.load_file('SM86-S-1.prmtop', 'SM86-S-1.pdb')

# Save a CHARMM PSF and crd file
amber.save('SM86-S-1.psf')
