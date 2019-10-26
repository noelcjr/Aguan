#!\bin\bash

a="/home/noel/Projects/FixMD/Aguan/build/jar/Aguan.jar"
echo "Finished Running tip3p t4"
echo "Start Running Reference minimization Tests."
echo "All other test will derive from this ones."
echo "Input file are the simme for thre types of waters,"
echo "I reuse them to save space and for consistency."
echo "density 1gr/cm^3"
echo "Directory Tests:"
echo "Electric_Fields  Temperature_Gradient  Thermostats  XYZ_Walls  Z_Walls"
cd Thermostats
java -jar $a -c -random tip3p_T1 1.8662 1.8662 1.8662 0 1 216 tip3p
java -jar $a -m tip3p_T1.in 0 500
java -jar $a -c -random tip4p_T1 1.8662 1.8662 1.8662 0 1 216 tip4p
java -jar $a -m tip4p_T1.in 0 500
java -jar $a -c -random tip5p_T1 1.8662 1.8662 1.8662 0 1 216 tip5p
java -jar $a -m tip5p_T1.in 0 500
java -jar $a -c -psf tip3p_T1_0.rst
mv tip3p_T1_0.psf ../
java -jar $a -c -psf tip4p_T1_0.rst
mv tip4p_T1_0.psf ../
java -jar $a -c -psf tip5p_T1_0.rst
mv tip5p_T1_0.psf ../
echo "End of initial minimizations"
echo "1  ######################################################"
java -jar $a -s tip3p_T1.in 1 100
java -jar $a -s tip4p_T1.in 1 100
java -jar $a -s tip5p_T1.in 1 100
echo "End thermostat 1 test"
echo "2  ######################################################"
cp tip3p_T1_1.rst tip3p_T2_1.rst
cp tip4p_T1_1.rst tip4p_T2_1.rst
cp tip5p_T1_1.rst tip5p_T2_1.rst
java -jar $a -s tip3p_T2.in 1 100
java -jar $a -s tip4p_T2.in 1 100
java -jar $a -s tip5p_T2.in 1 100
echo "End thermostat 2 test"
echo "3  ######################################################"
cp tip3p_T1_1.rst tip3p_T3_1.rst
cp tip4p_T1_1.rst tip4p_T3_1.rst
cp tip5p_T1_1.rst tip5p_T3_1.rst
java -jar $a -s tip3p_T3.in 1 100
java -jar $a -s tip4p_T3.in 1 100
java -jar $a -s tip5p_T3.in 1 100
echo "End thermostat 3 test"
echo "4 #######################################################"
cd ../Electric_Fields
cp ../Thermostats/tip3p_T1_1.rst tip3p_150EE_1.rst
cp ../Thermostats/tip4p_T1_1.rst tip4p_150EE_1.rst
cp ../Thermostats/tip5p_T1_1.rst tip5p_150EE_1.rst
java -jar $a -s tip3p_150EE.in 1 100
java -jar $a -s tip4p_150EE.in 1 100
java -jar $a -s tip5p_150EE.in 1 100
echo "End of Electric Field test"
echo "5 ########################################################"
cd ../Temperature_Gradient
cp ../Thermostats/tip3p_T1_1.rst tip3p_0.01TG_T2_1.rst
cp ../Thermostats/tip4p_T1_1.rst tip4p_0.01TG_T2_1.rst
cp ../Thermostats/tip5p_T1_1.rst tip5p_0.01TG_T2_1.rst
java -jar $a -s tip3p_0.01TG_T2.in 1 100
java -jar $a -s tip4p_0.01TG_T2.in 1 100
java -jar $a -s tip5p_0.01TG_T2.in 1 100
echo "End of Temperature Gradient test"
echo "6 ########################################################"
cd ../XYZ_Walls
java -jar $a -c -random tip3p_xyzWBC_0W 1.8662 1.8662 1.8662 0 1 216 tip3p
java -jar $a -m tip3p_xyzWBC_0W.in 0 500
java -jar $a -c -random tip4p_xyzWBC_0W 1.8662 1.8662 1.8662 0 1 216 tip4p
java -jar $a -m tip4p_xyzWBC_0W.in 0 500
java -jar $a -c -random tip5p_xyzWBC_0W 1.8662 1.8662 1.8662 0 1 216 tip5p
java -jar $a -m tip5p_xyzWBC_0W.in 0 500

java -jar $a -s tip3p_xyzWBC_0W.in 1 100
java -jar $a -s tip4p_xyzWBC_0W.in 1 100
java -jar $a -s tip5p_xyzWBC_0W.in 1 100
echo "End of xyz wall 0 Boundary Condition test"
echo "7  ######################################################"
cd ../Z_Walls
java -jar $a -c -random tip3p_zWBC_0W 1.8662 1.8662 1.8662 0 1 216 tip3p
java -jar $a -m tip3p_zWBC_0W.in 0 500
java -jar $a -c -random tip4p_zWBC_0W 1.8662 1.8662 1.8662 0 1 216 tip4p
java -jar $a -m tip4p_zWBC_0W.in 0 500
java -jar $a -c -random tip5p_zWBC_0W 1.8662 1.8662 1.8662 0 1 216 tip5p
java -jar $a -m tip5p_zWBC_0W.in 0 500

java -jar $a -s tip3p_zWBC_0W.in 1 100
java -jar $a -s tip4p_zWBC_0W.in 1 100
java -jar $a -s tip5p_zWBC_0W.in 1 100
echo "End of z wall 0 Boundary Condition test"
echo "8  ######################################################"
