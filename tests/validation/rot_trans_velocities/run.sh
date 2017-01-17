#!/bin/bash

a="/home/noel/Projects/FixMD/Aguan/build/jar/Aguan.jar"
echo "Generate RST File"
# Do not uncomment and run the following line. The restart file was modified manually for this test.
#java -jar $a -c -random tip3p_8 8 tip3p 2 2 2 0 > ang_vels.out
java -jar $a -c -psf tip3p_8_0.rst
java -jar $a -c -crd tip3p_8_0.rst
java -jar $a -s tip3p_8.in 0 3600
