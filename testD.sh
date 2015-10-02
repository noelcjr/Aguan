bash ant.sh
perl run.prl 1
mv out.txt out.txt.D
mv edr.out edr.out.D
mv trr.out trr.out.D
vimdiff edr.out ../Gromacs/tip3p/two_molecules_2ps/tip3p_run1.edr.out
