Install LOOS from https://github.com/GrossfieldLab/loos
Run the following command from terminal to get the cross-angles for the neutral ionizable lipids
cross-dist --num_bins 25 --cutoff 10 -s 'resname=="SN" || resname=="RN"' --use-cosine=1 -i 1000 -k 200000 parameter-file.prmtop trajectory.dcd
Run the following commant to get the cross-angles for the protonated ionizable lipids
cross-dist --num_bins 25 --cutoff 10 -s 'resname=="SP" || resname=="RP"' --use-cosine=1 -i 1000 -k 200000 parameter-file.prmtop trajectory.dcd
