package Aguan;

import Aguan.TheMatrix.TheMatrix;
import java.io.*;
import java.text.*;
/**
 * Time Step 0.00060241
 * @author Noel Carrascal
 */
public class parameters{
    public static int indexCount;
    private parameters(){
            indexCount = 0;
    }
    public static void DefinePar(TheMatrix TM, String name){
        TM.molType = name;
        if(name.equalsIgnoreCase("tip3p")){
           TM.sitesMol[TM.nMolTypes] = 4;
           TM.bCon[TM.nMolTypes] = 120.4926;             // 120.493009; //120.492627; // old 120.3995; 
           TM.ro[TM.nMolTypes] = 0.315061;               //  3.15061;     
           TM.ep[TM.nMolTypes] = 0.6364;                 // 0.6364 kj/mol x 1 kcal / 4.184 kj =  0.1521  kcal/mol; 
           TM.pre_convFact[TM.nMolTypes] = 33.789830;    // MegaPascals
           TM.density_convFact[TM.nMolTypes] = 1.043226; // for density 997047.9 g/m^3 at 25 Celcius density varies with temperature.
           TM.mInertX[TM.nMolTypes] = 0.00980;
           TM.mInertY[TM.nMolTypes] = 0.00340;
           TM.mInertZ[TM.nMolTypes] = 0.00640;        
        }else if(name.equalsIgnoreCase("tip4p")){
           TM.sitesMol[TM.nMolTypes] = 4;
           TM.bCon[TM.nMolTypes] = 183.8363;
           TM.ro[TM.nMolTypes] = 0.315365;
           TM.ep[TM.nMolTypes] = 0.6480;
           TM.pre_convFact[TM.nMolTypes] = 34.3070933;
           TM.density_convFact[TM.nMolTypes] = 1.046249;
           TM.mInertX[TM.nMolTypes] = 0.00980;
           TM.mInertY[TM.nMolTypes] = 0.00340;
           TM.mInertZ[TM.nMolTypes] = 0.00640;
        }else if(name.equalsIgnoreCase("tip5p")){
           TM.sitesMol[TM.nMolTypes] = 5;
           TM.bCon[TM.nMolTypes] = 38.63732;
           TM.ro[TM.nMolTypes] = 0.31200;
           TM.ep[TM.nMolTypes] = 0.6694;
           TM.pre_convFact[TM.nMolTypes] = 36.599177;
           TM.density_convFact[TM.nMolTypes] = 1.013114;
           TM.mInertX[TM.nMolTypes] = 0.00980;
           TM.mInertY[TM.nMolTypes] = 0.00340;
           TM.mInertZ[TM.nMolTypes] = 0.00640;
        }else if(name.equalsIgnoreCase("wt0LJ") || name.equalsIgnoreCase("wt1LJ")  || name.equalsIgnoreCase("wt2LJ")  || name.equalsIgnoreCase("wt3LJ")  || name.equalsIgnoreCase("wt4LJ")  || name.equalsIgnoreCase("wt5LJ")){
           // Modeled as a monoatomic tip3p with LJ and charge at the center
           TM.sitesMol[TM.nMolTypes] = 2;
           // If the waters are generated before the LJ particles then ro and bcon from the water
           // will be used in the calculations as a way to not have to modify the code anymore and be proned to errors.
           // This is a heack and easier way to do this to get results quickly, it will be corrected in future versions.
           TM.bCon[TM.nMolTypes] = 120.4926;             // 120.493009; //120.492627; // old 120.3995; 
           TM.ro[TM.nMolTypes] = 0.315061;               //  3.15061;     
           TM.ep[TM.nMolTypes] = 0.6364;                 // 0.6364 kj/mol x 1 kcal / 4.184 kj =  0.1521  kcal/mol; 
           TM.pre_convFact[TM.nMolTypes] = 33.789830;    // MegaPascals
           TM.density_convFact[TM.nMolTypes] = 1.043226; // for density 997047.9 g/m^3 at 25 Celcius density varies with temperature.
           TM.mInertX[TM.nMolTypes] = 0.0;
           TM.mInertY[TM.nMolTypes] = 0.0;
           TM.mInertZ[TM.nMolTypes] = 0.0;
           // Atomic gas with Lenard jones parameters identical to water.
           // This molecule is ficticious but it makes the modification of the code easier.
           // It is not much different from a noble gas in its parameters.
           // Noble gases will be added later when the code is made capable of dealing with
           // many atom types and their paramters.
        }else if(name.equalsIgnoreCase("wt6LJ") || name.equalsIgnoreCase("wt7LJ")  || name.equalsIgnoreCase("wt8LJ")  || name.equalsIgnoreCase("wt9LJ")  || name.equalsIgnoreCase("wtXLJ")){
           // Modeled as a monoatomic tip3p with LJ and charge at the center
           TM.sitesMol[TM.nMolTypes] = 2;
           // If the waters are generated before the LJ particles then ro and bcon from the water
           // will be used in the calculations as a way to not have to modify the code anymore and be proned to errors.
           // This is a heack and easier way to do this to get results quickly, it will be corrected in future versions.
           TM.bCon[TM.nMolTypes] = 120.4926;             // 120.493009; //120.492627; // old 120.3995; 
           TM.ro[TM.nMolTypes] = 0.315061;               //  3.15061;     
           TM.ep[TM.nMolTypes] = 0.6364;                 // 0.6364 kj/mol x 1 kcal / 4.184 kj =  0.1521  kcal/mol; 
           TM.pre_convFact[TM.nMolTypes] = 33.789830;    // MegaPascals
           TM.density_convFact[TM.nMolTypes] = 1.043226; // for density 997047.9 g/m^3 at 25 Celcius density varies with temperature.
           TM.mInertX[TM.nMolTypes] = 0.0;
           TM.mInertY[TM.nMolTypes] = 0.0;
           TM.mInertZ[TM.nMolTypes] = 0.0;
           // Atomic gas with Lenard jones parameters identical to water.
           // This molecule is ficticious but it makes the modification of the code easier.
           // It is not much different from a noble gas in its parameters.
           // Noble gases will be added later when the code is made capable of dealing with
           // many atom types and their paramters.
        }else if(name.equalsIgnoreCase("st2")){
           TM.sitesMol[TM.nMolTypes] = 5;
           TM.bCon[TM.nMolTypes] = 83.892285; // Fixed
           TM.ro[TM.nMolTypes] = 0.31000;  // Fixed
           TM.ep[TM.nMolTypes] = 0.31694;  // Fixed
           TM.pre_convFact[TM.nMolTypes] = 17.66611; // fixed
           TM.density_convFact[TM.nMolTypes] = 0.993755; // Fixed
           TM.mInertX[TM.nMolTypes] = 0.00980; // 
           TM.mInertY[TM.nMolTypes] = 0.00340; //
           TM.mInertZ[TM.nMolTypes] = 0.00640; //
        }else{}/*else if(name.equalsIgnoreCase("spc")){
        
        }else if(name.equalsIgnoreCase("ne")){
        */
    }
    public static void DefineMol(TheMatrix TM, String name){
        System.out.println("DefineMol ---- dummy"); 
    }
    public static void DefineMol(TheMatrix TM, int numbMolecs, String name){
        //for(int j = 0; j < TM.sitesMol[nMolNumber]; j++){
        //    TM.rmx[j] = 0.0; TM.rmy[j] = 0.0; TM.rmz[j] = 0.0;
        //}
        if(name.equalsIgnoreCase("tip3p")){
           for(int i = 0; i < numbMolecs; i++){
               TM.rmz[indexCount+0] =-0.0207;
               TM.rmz[indexCount+1] =-0.0207;
               TM.rmy[indexCount+2] = 0.2402;   TM.rmz[indexCount+2] = 0.1653;
               TM.rmy[indexCount+3] =-0.2402;   TM.rmz[indexCount+3] = 0.1653;
               TM.typeF[indexCount+0] = 1;      TM.typeF[indexCount+1] = 2;
               TM.typeF[indexCount+2] = 3;      TM.typeF[indexCount+3] = 3;
               TM.atomType[indexCount] = 1;
               TM.atomType[indexCount+1] = 2;
               TM.atomType[indexCount+2] = 3;
               TM.atomType[indexCount+3] = 3;
               TM.moleculeIndex[i] = indexCount;
               TM.moleculeIndex[i+1] = indexCount + TM.sitesMol[TM.nMolTypes];
               indexCount = indexCount + TM.sitesMol[TM.nMolTypes];
           }
           TM.nAtoms[TM.nMolTypes] = 3;     TM.nCharges[TM.nMolTypes] = 3;
           TM.nPoints[TM.nMolTypes] = 3;     TM.vdwPoint[TM.nMolTypes] = 1;
        }else if(name.equalsIgnoreCase("wt0LJ")){
           for(int i = 0; i < numbMolecs; i++){
               // TM.rmz[indexCount+0] = -1; 
               TM.typeF[indexCount+0] = 1;   TM.typeF[indexCount+1] = 10;
               TM.LJcharge = 0.0;
               TM.atomType[indexCount] = 1;
               TM.atomType[indexCount+1] = 2;
               TM.moleculeIndex[i] = indexCount;
               TM.moleculeIndex[i+1] = indexCount + TM.sitesMol[TM.nMolTypes]; 
               indexCount = indexCount + TM.sitesMol[TM.nMolTypes];
           }
           TM.nAtoms[TM.nMolTypes] = 1;     TM.nCharges[TM.nMolTypes] = 1;
           TM.nPoints[TM.nMolTypes] = 1;    TM.vdwPoint[TM.nMolTypes] = 1;
        }else if(name.equalsIgnoreCase("wt1LJ")){
           for(int i = 0; i < numbMolecs; i++){
               // TM.rmz[indexCount+0] = -1; 
               TM.typeF[indexCount+0] = 1;   TM.typeF[indexCount+1] = 10;
               TM.LJcharge = 0.1;
               TM.atomType[indexCount] = 1;
               TM.atomType[indexCount+1] = 2;
               TM.moleculeIndex[i] = indexCount;
               TM.moleculeIndex[i+1] = indexCount + TM.sitesMol[TM.nMolTypes];
               indexCount = indexCount + TM.sitesMol[TM.nMolTypes];
           }
           TM.nAtoms[TM.nMolTypes] = 1;     TM.nCharges[TM.nMolTypes] = 1;
           TM.nPoints[TM.nMolTypes] = 1;    TM.vdwPoint[TM.nMolTypes] = 1;
        }else if(name.equalsIgnoreCase("wt2LJ")){
           for(int i = 0; i < numbMolecs; i++){
               // TM.rmz[indexCount+0] = -1; 
               TM.typeF[indexCount+0] = 1;   TM.typeF[indexCount+1] = 10;
               TM.LJcharge = 0.2;
               TM.atomType[indexCount] = 1;
               TM.atomType[indexCount+1] = 2;
               TM.moleculeIndex[i] = indexCount;
               TM.moleculeIndex[i+1] = indexCount + TM.sitesMol[TM.nMolTypes];
               indexCount = indexCount + TM.sitesMol[TM.nMolTypes];
           }
           TM.nAtoms[TM.nMolTypes] = 1;     TM.nCharges[TM.nMolTypes] = 1;
           TM.nPoints[TM.nMolTypes] = 1;    TM.vdwPoint[TM.nMolTypes] = 1;
        }else if(name.equalsIgnoreCase("wt3LJ")){
           for(int i = 0; i < numbMolecs; i++){
               // TM.rmz[indexCount+0] = -1; 
               TM.typeF[indexCount+0] = 1;   TM.typeF[indexCount+1] = 10;
               TM.LJcharge = 0.3;
               TM.atomType[indexCount] = 1;
               TM.atomType[indexCount+1] = 2;
               TM.moleculeIndex[i] = indexCount;
               TM.moleculeIndex[i+1] = indexCount + TM.sitesMol[TM.nMolTypes];
               indexCount = indexCount + TM.sitesMol[TM.nMolTypes];
           }
           TM.nAtoms[TM.nMolTypes] = 1;     TM.nCharges[TM.nMolTypes] = 1;
           TM.nPoints[TM.nMolTypes] = 1;    TM.vdwPoint[TM.nMolTypes] = 1;
        }else if(name.equalsIgnoreCase("wt4LJ")){
           for(int i = 0; i < numbMolecs; i++){
               // TM.rmz[indexCount+0] = -1; 
               TM.typeF[indexCount+0] = 1;   TM.typeF[indexCount+1] = 10;
               TM.LJcharge = 0.0;
               TM.atomType[indexCount] = 1;
               TM.atomType[indexCount+1] = 2;
               TM.moleculeIndex[i] = indexCount;
               TM.moleculeIndex[i+1] = indexCount + TM.sitesMol[TM.nMolTypes];
               indexCount = indexCount + TM.sitesMol[TM.nMolTypes];
           }
           TM.nAtoms[TM.nMolTypes] = 1;     TM.nCharges[TM.nMolTypes] = 1;
           TM.nPoints[TM.nMolTypes] = 1;    TM.vdwPoint[TM.nMolTypes] = 1;
        }else if(name.equalsIgnoreCase("wt5LJ")){
           for(int i = 0; i < numbMolecs; i++){
               // TM.rmz[indexCount+0] = -1; 
               TM.typeF[indexCount+0] = 1;   TM.typeF[indexCount+1] = 10;
               TM.LJcharge = 0.5;
               TM.atomType[indexCount] = 1;
               TM.atomType[indexCount+1] = 2;
               TM.moleculeIndex[i] = indexCount;
               TM.moleculeIndex[i+1] = indexCount + TM.sitesMol[TM.nMolTypes];
               indexCount = indexCount + TM.sitesMol[TM.nMolTypes];
           }
           TM.nAtoms[TM.nMolTypes] = 1;     TM.nCharges[TM.nMolTypes] = 1;
           TM.nPoints[TM.nMolTypes] = 1;    TM.vdwPoint[TM.nMolTypes] = 1;
        }else if(name.equalsIgnoreCase("wt6LJ")){
           for(int i = 0; i < numbMolecs; i++){
               // TM.rmz[indexCount+0] = -1; 
               TM.typeF[indexCount+0] = 1;   TM.typeF[indexCount+1] = 10;
               TM.LJcharge = 0.6;
               TM.atomType[indexCount] = 1;
               TM.atomType[indexCount+1] = 2;
               TM.moleculeIndex[i] = indexCount;
               TM.moleculeIndex[i+1] = indexCount + TM.sitesMol[TM.nMolTypes];
               indexCount = indexCount + TM.sitesMol[TM.nMolTypes];
           }
           TM.nAtoms[TM.nMolTypes] = 1;     TM.nCharges[TM.nMolTypes] = 1;
           TM.nPoints[TM.nMolTypes] = 1;    TM.vdwPoint[TM.nMolTypes] = 1;
        }else if(name.equalsIgnoreCase("wt7LJ")){
           for(int i = 0; i < numbMolecs; i++){
               // TM.rmz[indexCount+0] = -1; 
               TM.typeF[indexCount+0] = 1;   TM.typeF[indexCount+1] = 10;
               TM.LJcharge = 0.7;
               TM.atomType[indexCount] = 1;
               TM.atomType[indexCount+1] = 2;
               TM.moleculeIndex[i] = indexCount;
               TM.moleculeIndex[i+1] = indexCount + TM.sitesMol[TM.nMolTypes];
               indexCount = indexCount + TM.sitesMol[TM.nMolTypes];
           }
           TM.nAtoms[TM.nMolTypes] = 1;     TM.nCharges[TM.nMolTypes] = 1;
           TM.nPoints[TM.nMolTypes] = 1;    TM.vdwPoint[TM.nMolTypes] = 1;
        }else if(name.equalsIgnoreCase("wt8LJ")){
           for(int i = 0; i < numbMolecs; i++){
               // TM.rmz[indexCount+0] = -1; 
               TM.typeF[indexCount+0] = 1;   TM.typeF[indexCount+1] = 10;
               TM.LJcharge = 0.8;
               TM.atomType[indexCount] = 1;
               TM.atomType[indexCount+1] = 2;
               TM.moleculeIndex[i] = indexCount;
               TM.moleculeIndex[i+1] = indexCount + TM.sitesMol[TM.nMolTypes];
               indexCount = indexCount + TM.sitesMol[TM.nMolTypes];
           }
           TM.nAtoms[TM.nMolTypes] = 1;     TM.nCharges[TM.nMolTypes] = 1;
           TM.nPoints[TM.nMolTypes] = 1;    TM.vdwPoint[TM.nMolTypes] = 1;
        }else if(name.equalsIgnoreCase("wt9LJ")){
           for(int i = 0; i < numbMolecs; i++){
               // TM.rmz[indexCount+0] = -1; 
               TM.typeF[indexCount+0] = 1;   TM.typeF[indexCount+1] = 10;
               TM.LJcharge = 0.9;
               TM.atomType[indexCount] = 1;
               TM.atomType[indexCount+1] = 2;
               TM.moleculeIndex[i] = indexCount;
               TM.moleculeIndex[i+1] = indexCount + TM.sitesMol[TM.nMolTypes];
               indexCount = indexCount + TM.sitesMol[TM.nMolTypes];
           }
           TM.nAtoms[TM.nMolTypes] = 1;     TM.nCharges[TM.nMolTypes] = 1;
           TM.nPoints[TM.nMolTypes] = 1;    TM.vdwPoint[TM.nMolTypes] = 1;
        }else if(name.equalsIgnoreCase("wtXLJ")){
           for(int i = 0; i < numbMolecs; i++){
               // TM.rmz[indexCount+0] = -1; 
               TM.typeF[indexCount+0] = 1;   TM.typeF[indexCount+1] = 10;
               TM.LJcharge = 1.0;
               TM.atomType[indexCount] = 1;
               TM.atomType[indexCount+1] = 2;
               TM.moleculeIndex[i] = indexCount;
               TM.moleculeIndex[i+1] = indexCount + TM.sitesMol[TM.nMolTypes];
               indexCount = indexCount + TM.sitesMol[TM.nMolTypes];
           }
           TM.nAtoms[TM.nMolTypes] = 1;     TM.nCharges[TM.nMolTypes] = 1;
           TM.nPoints[TM.nMolTypes] = 1;    TM.vdwPoint[TM.nMolTypes] = 1;
        }else if(name.equalsIgnoreCase("tip4p")){
           for(int i = 0; i < numbMolecs; i++){
               TM.rmz[indexCount+0] =-0.02065;
               TM.rmz[indexCount+1] = 0.02695;
               TM.rmy[indexCount+2] = 0.2400;    TM.rmz[indexCount+2] = 0.1652;
               TM.rmy[indexCount+3] = -0.2400;   TM.rmz[indexCount+3] = 0.1652;
               TM.typeF[indexCount+0] = 1;       TM.typeF[indexCount+1] = 2;
               TM.typeF[indexCount+2] = 3;       TM.typeF[indexCount+3] = 3;
               TM.atomType[indexCount] = 1;
               TM.atomType[indexCount+1] = 2;
               TM.atomType[indexCount+2] = 3;
               TM.atomType[indexCount+3] = 3;
               TM.moleculeIndex[i] = indexCount;
               TM.moleculeIndex[i+1] = indexCount + TM.sitesMol[TM.nMolTypes];
               indexCount = indexCount + TM.sitesMol[TM.nMolTypes];
           }
           TM.nAtoms[TM.nMolTypes] = 3;      TM.nCharges[TM.nMolTypes] = 3;
           TM.nPoints[TM.nMolTypes] = 3;     TM.vdwPoint[TM.nMolTypes] = 1;
        }else if(name.equalsIgnoreCase("tip5p")){
           for(int i = 0; i < numbMolecs; i++){
              TM.rmz[indexCount+0] = -0.0209;
              TM.rmy[indexCount+1] = 0.2426;   TM.rmz[indexCount+1] = 0.1669;
              TM.rmy[indexCount+2] =-0.2426;   TM.rmz[indexCount+2] = 0.1669;
              TM.rmx[indexCount+3] = 0.1832;   TM.rmz[indexCount+3] = -0.1504;
              TM.rmx[indexCount+4] =-0.1832;   TM.rmz[indexCount+4] = -0.1504;
              TM.typeF[indexCount+0] = 1;      TM.typeF[indexCount+1] = 3;    
              TM.typeF[indexCount+2] = 3;      TM.typeF[indexCount+3] = 7;       TM.typeF[indexCount+4] = 7;
              TM.atomType[indexCount] = 2;
              TM.atomType[indexCount+1] = 3;
              TM.atomType[indexCount+2] = 3;
              TM.atomType[indexCount+3] = 4;
              TM.atomType[indexCount+4] = 4;
              TM.moleculeIndex[i] = indexCount;
              TM.moleculeIndex[i+1] = indexCount + TM.sitesMol[TM.nMolTypes];
              indexCount = indexCount + TM.sitesMol[TM.nMolTypes];
           }
           TM.nAtoms[TM.nMolTypes] = 3;     TM.nCharges[TM.nMolTypes] = 4;
           TM.nPoints[TM.nMolTypes] = 5;    TM.vdwPoint[TM.nMolTypes] = 7;
        }else if(name.equalsIgnoreCase("st2")){
           for(int i = 0; i < numbMolecs; i++){
              TM.rmx[indexCount+0] = -0.0207;
              TM.rmy[indexCount+1] = 0.263;   TM.rmz[indexCount+1] = 0.1656;
              TM.rmy[indexCount+2] =-0.263;   TM.rmz[indexCount+2] = 0.1656;
              TM.rmx[indexCount+3] = 0.2107;  TM.rmz[indexCount+3] =-0.1697;
              TM.rmx[indexCount+4] =-0.2107;  TM.rmz[indexCount+4] =-0.1697;
              TM.typeF[indexCount+0] = 1;     TM.typeF[indexCount+1] = 3;    
              TM.typeF[indexCount+2] = 3;     TM.typeF[indexCount+3] = 7;   TM.typeF[indexCount+4] = 7;
              TM.atomType[indexCount] = 2;
              TM.atomType[indexCount+1] = 3;
              TM.atomType[indexCount+2] = 3;
              TM.atomType[indexCount+3] = 4;
              TM.atomType[indexCount+4] = 4;
              TM.moleculeIndex[i] = indexCount; 
              TM.moleculeIndex[i+1] = indexCount + TM.sitesMol[TM.nMolTypes];
              indexCount = indexCount + TM.sitesMol[TM.nMolTypes];
           }
           TM.nAtoms[TM.nMolTypes] = 3;     TM.nCharges[TM.nMolTypes] = 4;
           TM.nPoints[TM.nMolTypes] = 5;    TM.vdwPoint[TM.nMolTypes] = 7;
        }else if(name.equalsIgnoreCase("spc")){

        }
    }
}
