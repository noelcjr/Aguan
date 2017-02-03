/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package Aguan.TheMatrix;

import Aguan.files.input;
import Aguan.files.restart;
import Aguan.parameters;
/**
 *
 * @author Noel
 */
public class TheMatrix {
    /////  Parameters from input file /////
    public double deltaT;                //
    public String boundaryCondition;     //
    public double rCut;                  //
    public int stepAvg;                  //
    public int stepLimit;                //
    public boolean moreCycles;           //
    public int trajOut;                  //
    public double targetTemperature;     //
    public double tau;                   //
    public double electricField;         //
    public double TempRateChange;        //
    public int ThermostatType;           // 
    public int wallType;                 //
    ///////////////////////////////////////
    /// Parameters from rst file header  ////////////
    public int stepCount;                          //
    public int restartCount;                       //
    public double density;                         //
    public int initUcellX, initUcellY, initUcellZ; //
    public double regionX, regionY, regionZ;       //
    public double volume;                          //
    public int nMol, nMolTypes;                    //
    public String[] nMolNames;
    public int totalSitesNumber, sitesNoVDW;       //
    ///////////////////////////////////////////////// 
    //// parameters in parameters.definePar /////
    public int[] sitesMol, sitesMolIdx;        //
    public double[] bCon;                      //  Conv Factor
    public double[] ro;                        //  Conv Factor
    public double[] ep;                        //  Conv Factor
    public double[] pre_convFact;              //  Conv Factor
    public double[] density_convFact;          //  Conv Factor
    public double[] mInertX, mInertY, mInertZ; //
    /////////////////////////////////////////////
    ///// Parameters in TM.initMatrix //////////////////////////////
    public int NDIM;                                              //
    public double velMag;                                         //
    public int[] atomType, moleculeIndex, moleculeTypeIndex;      //
    public double[] rMatT;                                        //
    public double[] rxs, rys, rzs, vxs, vys, vzs, fxs, fys, fzs;  //
    public double[] rx, ry, rz, rvx, rvy, rvz;                    //
    public double[] rax, ray, raz, wax, way, waz, wvx, wvy, wvz;  //
    public double[] rmx, rmy, rmz, vmx, vmy, vmz;                 //
    public double[] q_u1, q_u2, q_u3, q_u4;                       //
    public double[] qv_u1, qv_u2, qv_u3, qv_u4;                   //
    public double[] av_1, av_2, av_3;                             //
    public int[] typeF;                                           //
    public int q, q0, q1, q2, q3, q4, q5, q10, q21, q30, q32, q40;//
    public int q41, q42, q43, q50, q51, q52, q53, q410, q421;     //
    public int q430, q432, q510, q521, q530, q532;                //
    ////////////////////////////////////////////////////////////////
    public boolean ReadORGen, changeDensity;
    public double oldDensity, pressure, pressureSum;
    public double instantTemperature, aveTemperature, rotationalTemperature, translationalTemperature;
    public double rotTempAve, trzTempAve, rotFrac, tranFrac;
    public int stepAdjustTemp, stepEquil;
    public String restartIn, restartOut, trajectory, psf, crd, logfile, outfile;
    public int[] nCharges, nAtoms, nPoints, vdwPoint, vdwPointIdx;
    public int allAtoms;
    public double vSumX, vSumY, vSumZ;
    public double wvSumX, wvSumY, wvSumZ;
    public double vvSum, vvrSum, virSum;
    public double uSum, uSumVDW, uSumEE, uSumRF1, uSumRF2;
    public double[] rvxf, rvyf, rvzf, rrvxf, rrvyf, rrvzf;
    public double totEnergyVal, totEnergySum, potEnergyVal, potEnergySum, kinEnergyVal, kinEnergySum;
    public double potEnergyVdwVal, potEnergyVdwSum, potEnergyEeVal, potEnergyEeSum;
    public double potEnergyRf1Val, potEnergyRf1Sum, potEnergyRf2Val, potEnergyRf2Sum, potEnergyEfVal, potEnergyEfSum;
    public double rotKinEnergyVal, rotKinEnergySum, trzKinEnergyVal, trzKinEnergySum;
    public String molType;
    public int minimization;
    /** Creates a new instance of TheMatrix */
    private input IN;
    private restart RR;
    public double LJcharge;
    public TheMatrix(){
       totalSitesNumber = 0;
       sitesNoVDW = 0;
       nMolTypes = 0;
       LJcharge = 1.0;
    }
    public void initParams(int numMolecules){
       sitesMol = new int[numMolecules];
       bCon = new double[numMolecules];
       ro = new double[numMolecules];
       ep = new double[numMolecules];
       pre_convFact = new double[numMolecules];
       density_convFact  = new double[numMolecules];
       mInertX = new double[numMolecules];
       mInertY = new double[numMolecules];
       mInertZ = new double[numMolecules];
       nAtoms = new int[numMolecules];
       nCharges = new int[numMolecules];
       nPoints = new int[numMolecules];
       vdwPoint = new int[numMolecules];
       moleculeTypeIndex = new int[numMolecules*2];
    }
    public void initMatrix(){
        NDIM = 3;
        //velMag = Math.sqrt(NDIM*(1.0-1.0/nMol)*targetTemperature);
        //replaced with this one to match previous version. But I am not sure it is ready.
        velMag = Math.sqrt(NDIM*(1.0-1.0/nMol)*targetTemperature*0.5); 
        atomType = new int[totalSitesNumber];
        sitesMolIdx = new int[nMol];
        vdwPointIdx = new int[nMol];
        moleculeIndex = new int[nMol*2];
        nMolNames = new String[nMol];
        rMatT = new double[9*nMol];

        typeF = new int[totalSitesNumber]; 
        rxs = new double[totalSitesNumber];  vxs = new double[totalSitesNumber];  fxs = new double[totalSitesNumber];
        rys = new double[totalSitesNumber];  vys = new double[totalSitesNumber];  fys = new double[totalSitesNumber];
        rzs = new double[totalSitesNumber];  vzs = new double[totalSitesNumber];  fzs = new double[totalSitesNumber];

        rmx = new double[totalSitesNumber];  vmx = new double[totalSitesNumber];
        rmy = new double[totalSitesNumber];  vmy = new double[totalSitesNumber];
        rmz = new double[totalSitesNumber];  vmz = new double[totalSitesNumber];
        for(int i = 0; i < totalSitesNumber; i++){
           rmx[i] = 0.0;  rmy[i] = 0.0; rmz[i] = 0.0;
        }
        rx =  new double[nMol];  rvx =  new double[nMol];
        ry =  new double[nMol];  rvy =  new double[nMol];
        rz =  new double[nMol];  rvz =  new double[nMol];

        rax = new double[nMol];  wax = new double[nMol];  wvx = new double[nMol];
        ray = new double[nMol];  way = new double[nMol];  wvy = new double[nMol];
        raz = new double[nMol];  waz = new double[nMol];  wvz = new double[nMol];

        q_u1 = new double[nMol];        qv_u1 = new double[nMol];    av_1 = new double[nMol];
        q_u2 = new double[nMol];        qv_u2 = new double[nMol];    av_2 = new double[nMol];
        q_u3 = new double[nMol];        qv_u3 = new double[nMol];    av_3 = new double[nMol];
        q_u4 = new double[nMol];        qv_u4 = new double[nMol];
        //typeF = new int[totalSitesNumber];
        q = q0 = q1 = q2 = q3 = q4 = q5 = q10 = q21 = q30 = q32 = q40 = q41 = q42 = q43 = q50 = q51 = q52 = q53 = q410 = q421 = q430 = q432 = q510 = q521 = q530 = q532 = 0;
    }
}
