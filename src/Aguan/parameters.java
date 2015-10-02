package Aguan;

import Aguan.TheMatrix.TheMatrix;
import java.io.*;
import java.text.*;
/**
 * Time Step 0.00060241
 * @author Noel Carrascal
 */
public class parameters{
    private parameters(){}
    public static void DefinePar(TheMatrix TM, String name){
        TM.molType = name;
        if(name.equalsIgnoreCase("tip3p")){
           TM.sitesMol = 4;
           TM.bCon = 120.4926;             // 120.493009; //120.492627; // old 120.3995; 
           TM.ro = 0.315061;               //  3.15061;     
           TM.ep = 0.6364;                 // 0.6364 kj/mol x 1 kcal / 4.184 kj =  0.1521  kcal/mol; 
           TM.pre_convFact = 33.789830;    // MegaPascals
           TM.density_convFact = 1.043226; // for density 997047.9 g/m^3 at 25 Celcius density varies with temperature.
           TM.mInertX = 0.00980;
           TM.mInertY = 0.00340;
           TM.mInertZ = 0.00640;        
        }else if(name.equalsIgnoreCase("tip4p")){
           TM.sitesMol = 4;
           TM.bCon = 183.8363;
           TM.ro = 0.315365;
           TM.ep = 0.6480;
           TM.pre_convFact = 34.3070933;
           TM.density_convFact = 1.046249;
           TM.mInertX = 0.00980;
           TM.mInertY = 0.00340;
           TM.mInertZ = 0.00640;
        }else if(name.equalsIgnoreCase("tip5p")){
           TM.sitesMol = 5;
           TM.bCon = 38.63732;
           TM.ro = 0.31200;
           TM.ep = 0.6694;
           TM.pre_convFact = 36.599177;
           TM.density_convFact = 1.013114;
           TM.mInertX = 0.00980;
           TM.mInertY = 0.00340;
           TM.mInertZ = 0.00640;
        }else{}/*else if(name.equalsIgnoreCase("st2")){
           TM.sitesMol = 5;
           TM.bCon = 83.8262;
           TM.ro = 3.1000;
           TM.ep =0.0757;
           TM.pre_convFact = 17.654342;
           TM.density_convFact = 0.9966979;
           TM.mInertX = 0.00980;
           TM.mInertY = 0.00340;
           TM.mInertZ = 0.00640;
        }else if(name.equalsIgnoreCase("spc")){
        
        }else if(name.equalsIgnoreCase("ne")){
         */
    }
    public static void DefineMol(TheMatrix TM, String name){
        for(int j = 0; j < TM.sitesMol; j++){
            TM.rmx[j] = 0.0; TM.rmy[j] = 0.0; TM.rmz[j] = 0.0;
        }
        if(name.equalsIgnoreCase("tip3p")){
           TM.rmz[0] = -0.0207;
           TM.rmz[1] = -0.0207;
           TM.rmy[2] = 0.2402;   TM.rmz[2] = 0.1653;
           TM.rmy[3] =-0.2402;   TM.rmz[3] = 0.1653;
           TM.typeF[0] = 1;   TM.typeF[1] = 2;    TM.typeF[2] = 3;    TM.typeF[3] = 3;
           TM.nAtoms = 3;     TM.nCharges = 3;    TM.nPoints = 3;     TM.vdwPoint = 1;
           // Only for tip4p it would have to be modified for other molecules or combination of molecules
           for(int j = 0; j < TM.nMol*TM.sitesMol; j=j+4){
               TM.atomType[j] = 1;
               TM.atomType[j+1] = 2;
               TM.atomType[j+2] = 3;
               TM.atomType[j+3] = 3;
           }
           for(int j = 0; j < TM.nMol; j++){TM.moleculeIndex[j] = (TM.sitesMol*3*j);} 
        }else if(name.equalsIgnoreCase("tip4p")){
           TM.rmz[0] =-0.02065;
           TM.rmz[1] = 0.02695;
           TM.rmy[2] = 0.2400;    TM.rmz[2] = 0.1652;
           TM.rmy[3] = -0.2400;   TM.rmz[3] = 0.1652;
           TM.typeF[0] = 1;   TM.typeF[1] = 2;    TM.typeF[2] = 3;    TM.typeF[3] = 3;
           TM.nAtoms = 3;      TM.nCharges = 3;   TM.nPoints = 3;     TM.vdwPoint = 1;
           // Only for tip4p it would have to be modified for other molecules or combination of molecules
           for(int j = 0; j < TM.nMol*TM.sitesMol; j=j+4){
               TM.atomType[j] = 1;
               TM.atomType[j+1] = 2;
               TM.atomType[j+2] = 3;
               TM.atomType[j+3] = 3;
           }
           for(int j = 0; j < TM.nMol; j++){TM.moleculeIndex[j] = (TM.sitesMol*3*j);}
        }else if(name.equalsIgnoreCase("tip5p")){
           TM.rmz[0] = -0.0209;
           TM.rmy[1] = 0.2426;   TM.rmz[1] = 0.1669;
           TM.rmy[2] =-0.2426;   TM.rmz[2] = 0.1669;
           TM.rmx[3] = 0.1832;   TM.rmz[3] = -0.1504;
           TM.rmx[4] =-0.1832;   TM.rmz[4] = -0.1504;
           TM.typeF[0] = 1;   TM.typeF[1] = 3;    TM.typeF[2] = 3;     TM.typeF[3] = 7;   TM.typeF[4] = 7;
           TM.nAtoms = 3;     TM.nCharges = 4;    TM.nPoints = 5;      TM.vdwPoint = 7;
           // Only for tip4p it would have to be modified for other molecules or combination of molecules
           for(int j = 0; j < TM.nMol*TM.sitesMol; j=j+5){
               TM.atomType[j] = 2;
               TM.atomType[j+1] = 3;
               TM.atomType[j+2] = 3;
               TM.atomType[j+3] = 4;
               TM.atomType[j+4] = 4;
           }
           for(int j = 0; j < TM.nMol; j++){TM.moleculeIndex[j] = (TM.sitesMol*3*j);}
        }else if(name.equalsIgnoreCase("st2")){
           TM.rmx[0] = -0.0207;
           TM.rmy[1] = 0.263;   TM.rmz[1] = 0.1656;
           TM.rmy[2] =-0.263;   TM.rmz[2] = 0.1656;
           TM.rmx[3] = 0.2107;  TM.rmz[3] =-0.1697;
           TM.rmx[4] =-0.2107;  TM.rmz[4] =-0.1697;
           TM.typeF[0] = 1;   TM.typeF[1] = 3;    TM.typeF[2] = 3;    TM.typeF[3] = 7;   TM.typeF[4] = 7;
           TM.nAtoms = 3;     TM.nCharges = 4;    TM.nPoints = 5;     TM.vdwPoint = 7;
           // Only for tip4p it would have to be modified for other molecules or combination of molecules
           for(int j = 0; j < TM.nMol*TM.sitesMol; j=j+5){
               TM.atomType[j] = 2;
               TM.atomType[j+1] = 3;
               TM.atomType[j+2] = 3;
               TM.atomType[j+3] = 4;
               TM.atomType[j+4] = 4;
           }
           for(int j = 0; j < TM.nMol; j++){TM.moleculeIndex[j] = (TM.sitesMol*3*j);}
        }else if(name.equalsIgnoreCase("spc")){

        }
    }
}
