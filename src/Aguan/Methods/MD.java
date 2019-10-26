package Aguan.Methods;

import Aguan.TheMatrix.TheMatrix;
import Aguan.files.restart;
import Aguan.IO.dcdManager;
import Aguan.parameters;
import java.io.*;
import java.text.*;
import Aguan.files.*;
import java.util.*;
/**
 * Time Step 0.00060241
 * @author Noel Carrascal
 */
public class MD{
    public log LOG;        public output OUT;
    public dcdManager dcdM;
    public TheMatrix TM;   public restart RR;    private input IN;
    double timeNow;        public String fileName;
    public MD(){}
    public MD(String[] args){
        TM = new TheMatrix();
        IN = new input(args[1]);
        IN.GetNameList(TM);
        fileName = args[1].substring(0,args[1].length()-3);
        TM.restartIn = fileName+"_"+args[2]+".rst";
        TM.stepLimit = Integer.parseInt(args[3]);
        RR = new restart(TM.restartIn);
        TM.targetTemperature = TM.targetTemperature*(3.8/298.0);///2;used2b div by 2
        RR.ReadRestartHeader(TM);
        TM.outfile = fileName+"_"+TM.restartCount+".out";       TM.trajectory = fileName+"_"+TM.restartCount+".dcd";
        TM.restartCount++;
        TM.restartOut = fileName+"_"+TM.restartCount+".rst";
        RR.setOutputFile(TM.restartOut);
        //TM.targetTemperature = TM.targetTemperature*(3.8/298.0);///2;used2b div by 2
        //TM.initMatrix(); 
        //parameters.DefineMol(TM,TM.molType);
        TM.deltaT = TM.deltaT*0.000625;
        // conversion to reduced units.//////////////////
        TM.rCut = TM.rCut/TM.ro[0];                    //
        /////////////////////////////////////////////////
        // TM.targetTemperature must be before TM.initMatrix
        // I want to move unit convertion to file/input.java
        //TM.targetTemperature = TM.targetTemperature*(3.8/298.0);///2;used2b div by 2
        TM.TempRateChange = TM.TempRateChange*(3.8/298.0);   
        TM.electricField = TM.electricField*(0.000007037/150);
        System.out.printf("nMols = %d ",TM.nMol);
        System.out.printf("boundaryCondition = %s, velMag = %10.5f\n",TM.boundaryCondition,TM.velMag);
        RR.ReadRestartFile(TM);
        OUT = new output(TM.outfile);
        dcdM = new dcdManager(TM.trajectory);
        GenSiteCoord(TM);
        if(args[0].equals("-m")){
           dcdM.write_header("CORD",((TM.stepLimit/TM.trajOut)+1),1, 1, TM.sitesNoVDW);
        }else if(args[0].equals("-s")){
           dcdM.write_header("CORD",(TM.stepLimit/TM.trajOut),1, 1, TM.sitesNoVDW);
        }
        if(TM.stepLimit > 0){TM.moreCycles = true;}else{TM.moreCycles = false;}
    }
    public void GenSiteCoord(TheMatrix TM){
        double tx, ty, tz;
        int index = 0;
        int atomCount = 1;
        int atomIndex = 0;
        int molCount = 1;
        for(int n = 0; n < TM.nMol; n++){
              for(int j = 0; j < TM.sitesMolIdx[n]; j++){
                tx = TM.rMatT[(n*9)+0]*TM.rmx[j] + TM.rMatT[(n*9)+3]*TM.rmy[j] + TM.rMatT[(n*9)+6]*TM.rmz[j];
                ty = TM.rMatT[(n*9)+1]*TM.rmx[j] + TM.rMatT[(n*9)+4]*TM.rmy[j] + TM.rMatT[(n*9)+7]*TM.rmz[j];
                tz = TM.rMatT[(n*9)+2]*TM.rmx[j] + TM.rMatT[(n*9)+5]*TM.rmy[j] + TM.rMatT[(n*9)+8]*TM.rmz[j];

                TM.rxs[index + j] = TM.rx[n] + tx;
                TM.rys[index + j] = TM.ry[n] + ty;
                TM.rzs[index + j] = TM.rz[n] + tz;
                atomCount++;
                atomIndex++;
              }
              index = index + TM.sitesMolIdx[n];
              molCount++;
        }
    }

    public void ComputeWallForcesXYZ(TheMatrix TM){
           int ms1F, ms1T;
           ms1F = ms1T = 0;
           double dr, rri, rri3, fcValx, fcValy, fcValz, wallCutOff, beginning_dist_4_wall_repulsion;
           // The effective size of the box is reduce by the wall potential because water does not 
           // even get close to the walls. A water should get as close to the wall as it is allowed
           // by the distance from a hydrogen nuclei to its effective radius. 
           // Hydrogen has a radius of 0.02245 nm in the charmm force field for TIP3P hydrogen.
           // I fixed this by adding TM.ro[0] to the location of the walls
           // WallCutOff considers molecules that are only in the lennard-jones repulsive limit.
           // This way, the walls effectively just acts only repealing waters and not attracting them.
           wallCutOff = 1.122462;//*TM.ro[0]; // = (2^(1/6))*TM.ro;
           beginning_dist_4_wall_repulsion = 0.74694;
           fcValx = fcValy = fcValz = 0;
           if(TM.wallType == 0){
              for(int m1 = 0; m1 < TM.nMol; m1++){
                  if((TM.regionX-TM.rx[m1]) < wallCutOff){ // Checks that you are within wall interaction rage.
                        dr = TM.regionX - TM.rx[m1] + beginning_dist_4_wall_repulsion;    // so that only repulsive forces are considered. 
                        rri = 1/(dr*dr);
                        rri3 = rri*rri*rri;
                        fcValx = -48 * rri3 * (rri3 - 0.5) *rri *dr;
                        TM.uSum += 4 * rri3 * (rri3 - 1) + 1;
                  }else if(TM.rx[m1] < wallCutOff){
                        dr = TM.rx[m1] + beginning_dist_4_wall_repulsion;
                        rri = 1/(dr*dr);
                        rri3 = rri*rri*rri;
                        fcValx = 48 * rri3 * (rri3 - 0.5) *rri *dr;
                        TM.uSum += 4 * rri3 * (rri3 - 1) + 1;
                  }else{// Need to deal with possibility of a molecule crossing the wall and be outside the box
                  }
                  if((TM.regionY-TM.ry[m1]) < wallCutOff){
                        dr = TM.regionY - TM.ry[m1] + beginning_dist_4_wall_repulsion;
                        rri = 1/(dr*dr);
                        rri3 = rri*rri*rri;
                        fcValy = -48 * rri3 * (rri3 - 0.5) *rri *dr;
                        TM.uSum += 4 * rri3 * (rri3 - 1) + 1;
                  }else if(TM.ry[m1] < wallCutOff){
                        dr = TM.ry[m1] + beginning_dist_4_wall_repulsion;
                        rri = 1/(dr*dr);
                        rri3 = rri*rri*rri;
                        fcValy = 48 * rri3 * (rri3 - 0.5) *rri *dr;
                        TM.uSum += 4 * rri3 * (rri3 - 1) + 1;
                  }else{// Need to deal with possibility of a molecule crossing the wall and be outside the box
                  }
                  if((TM.regionZ-TM.rz[m1]) < wallCutOff){
                        dr = TM.regionZ - TM.rz[m1] + beginning_dist_4_wall_repulsion;
                        rri = 1/(dr*dr);
                        rri3 = rri*rri*rri;
                        fcValz = -48 * rri3 * (rri3 - 0.5) *rri *dr;
                        TM.uSum += 4 * rri3 * (rri3 - 1) + 1;
                  }else if(TM.rz[m1] < wallCutOff){
                        dr = TM.rz[m1] + beginning_dist_4_wall_repulsion;
                        rri = 1/(dr*dr);
                        rri3 = rri*rri*rri;
                        fcValz = 48 * rri3 * (rri3 - 0.5) *rri *dr;
                        TM.uSum += 4 * rri3 * (rri3 - 1) + 1;
                  }else{// Need to deal with possibility of a molecule crossing the wall and be outside the box
                  }
                  ms1T = ms1T + TM.sitesMolIdx[m1];
                  //ms1 = m1 * TM.sitesMol[0];
                  for(int j1 = ms1F; j1 < ms1T; j1++){
                      TM.fxs[j1] = TM.fxs[j1] + fcValx;
                      TM.fys[j1] = TM.fys[j1] + fcValy;
                      TM.fzs[j1] = TM.fzs[j1] + fcValz;
                  }
                  TM.virSum = TM.virSum + (fcValx + fcValy + fcValz);
                  ms1F = ms1T;
              }
          }
    }
    public void ComputeWallForcesZ(TheMatrix TM){
           int ms1F, ms1T;
           ms1F = ms1T = 0;
           double dr, rri, rri3, fcValz, wallCutOff, beginning_dist_4_wall_repulsion;
           wallCutOff = 1.122462; //*TM.ro[0]; // = (2^(1/6))*TM.ro;
           beginning_dist_4_wall_repulsion = 0.74694;
           fcValz = 0;
           if(TM.wallType == 0){
              for(int m1 = 0; m1 < TM.nMol; m1++){
                  if((TM.regionZ-TM.rz[m1]) < wallCutOff){
                        dr = TM.regionZ - TM.rz[m1] + beginning_dist_4_wall_repulsion;
                        rri = 1/(dr*dr);
                        rri3 = rri*rri*rri;
                        fcValz = -48 * rri3 * (rri3 - 0.5) *rri *dr;
                        TM.uSum += 4 * rri3 * (rri3 - 1) + 1;
                  }else if(TM.rz[m1] < wallCutOff){
                        dr = TM.rz[m1] + beginning_dist_4_wall_repulsion;
                        rri = 1/(dr*dr);
                        rri3 = rri*rri*rri;
                        fcValz = 48 * rri3 * (rri3 - 0.5) *rri *dr;
                        TM.uSum += 4 * rri3 * (rri3 - 1) + 1;
                  }else{// Need to deal with possibility of a molecule crossing the wall and be outside the box
                  }
                  // Possible Error, or find out why I looped twice for the Z wall.
                  //ms1 = m1 * TM.sitesMol[0];
                  //for(int j1 = 0; j1 < TM.sitesMol[0] ; j1++){
                  //    TM.fzs[ms1+j1] = TM.fzs[ms1+j1] + fcValz;
                  //}
                  ms1T = ms1T + TM.sitesMolIdx[m1];
                  for(int j1 = ms1F; j1 < ms1T; j1++){
                      TM.fzs[j1] = TM.fzs[j1] + fcValz;
                  }
                  TM.virSum = TM.virSum + fcValz;
                  ms1F = ms1T;
              }
         }
    }

    public void ApplyBoundaryCond(TheMatrix TM){
    // Boundary conditions in the precense of walls may be different due to some particles
    // crossing the wall because they are not the center of mass used for image convertion.
    //   if(TM.boundaryCondition == 0){     // Withou boundary conditions and an electric field
                                            // a twister/free fall effect is observed.
    //   }else if(TM.boundaryCondition == 1){
         for(int n = 0; n < TM.nMol; n++){
                if( TM.rx[n] >= TM.regionX){TM.rx[n] = TM.rx[n] - TM.regionX;}
           else if( TM.rx[n] < 0){TM.rx[n] = TM.rx[n] + TM.regionX;}
                if( TM.ry[n] >= TM.regionY){TM.ry[n] = TM.ry[n] - TM.regionY;}
           else if( TM.ry[n] < 0){TM.ry[n] = TM.ry[n] + TM.regionY;}
                if( TM.rz[n] >= TM.regionZ){TM.rz[n] = TM.rz[n] - TM.regionZ;}
           else if( TM.rz[n] < 0){TM.rz[n] = TM.rz[n] + TM.regionZ;}
         }
     //  }
    }
    public void EvalProps(TheMatrix TM){
         double wvBx, wvBy, wvBz;
         wvBx = wvBy = wvBz = 0;
         TM.vSumX = TM.vSumY = TM.vSumZ = 0.0;
         TM.vvSum = TM.vvrSum = 0.0;
         double TKE = 0.0; double WKE = 0.0;
         for(int n = 0; n < TM.nMol; n++){
             TM.vSumX = TM.vSumX + TM.rvx[n];   // This comented lines could be used for
             TM.vSumY = TM.vSumY + TM.rvy[n];   // correcting drifting in MD and possibly
             TM.vSumZ = TM.vSumZ + TM.rvz[n];   // to avoid the 'flying ice problem'
             TM.vvSum += ((TM.rvx[n]*TM.rvx[n]) + (TM.rvy[n]*TM.rvy[n]) + (TM.rvz[n]*TM.rvz[n]));
             wvBx = TM.rMatT[(n*9)+0]*TM.wvx[n] + TM.rMatT[(n*9)+1]*TM.wvy[n] + TM.rMatT[(n*9)+2]*TM.wvz[n];
             wvBy = TM.rMatT[(n*9)+3]*TM.wvx[n] + TM.rMatT[(n*9)+4]*TM.wvy[n] + TM.rMatT[(n*9)+5]*TM.wvz[n];  
             wvBz = TM.rMatT[(n*9)+6]*TM.wvx[n] + TM.rMatT[(n*9)+7]*TM.wvy[n] + TM.rMatT[(n*9)+8]*TM.wvz[n];
             TM.vvrSum += ((TM.mInertX[0]*wvBx*wvBx) + (TM.mInertY[0]*wvBy*wvBy) + (TM.mInertZ[0]*wvBz*wvBz));                     
         }
         //System.out.println("EvalProps   TM.vvSum="+TM.vvSum+"    TM.vvrSum="+TM.vvrSum);
         TM.translationalTemperature = (TM.vvSum/TM.nMol)/3.0;
         TM.rotationalTemperature = (TM.vvrSum/TM.nMol)/3.0;
         TM.instantTemperature = TM.translationalTemperature + TM.rotationalTemperature;
         TM.rotKinEnergyVal = 0.5 * TM.vvrSum/TM.nMol;
         TM.trzKinEnergyVal = 0.5 * TM.vvSum/TM.nMol;
         TM.kinEnergyVal = TM.rotKinEnergyVal + TM.trzKinEnergyVal;
         TM.pressure = 2*((TM.vvSum + TM.vvrSum + TM.virSum)/TM.nMol)/(TM.NDIM*TM.volume);
         //Alternative way to measure P. TM.pressure = TM.density*((TM.vvSum + TM.vvrSum) + TM.virSum)/(TM.NDIM*TM.nMol);
         TM.potEnergyVdwVal = TM.uSumVDW/TM.nMol;
         TM.potEnergyEeVal = TM.uSumEE/TM.nMol;
         TM.potEnergyRf1Val = TM.uSumRF1/TM.nMol;
         TM.potEnergyRf2Val = TM.uSumRF2/TM.nMol;         
         TM.potEnergyVal = TM.potEnergyVdwVal + TM.potEnergyEeVal + TM.potEnergyRf1Val + TM.potEnergyRf2Val + TM.potEnergyEfVal;
         TM.totEnergyVal = TM.kinEnergyVal + TM.potEnergyVal;
         TM.rotFrac = 0.5;
         TM.tranFrac = 0.5;
         if(TM.ThermostatType == 4){
             TM.rotFrac = TM.rotationalTemperature/TM.instantTemperature;
             TM.tranFrac = TM.translationalTemperature/TM.instantTemperature;
             for(int n = 0; n < TM.nMol; n++){
                 TM.rvxf[n] = TM.rvx[n]/TM.vvSum;
                 TM.rvyf[n] = TM.rvy[n]/TM.vvSum;
                 TM.rvzf[n] = TM.rvz[n]/TM.vvSum;
                 wvBx = TM.rMatT[(n*9)+0]*TM.wvx[n] + TM.rMatT[(n*9)+1]*TM.wvy[n] + TM.rMatT[(n*9)+2]*TM.wvz[n];
                 wvBy = TM.rMatT[(n*9)+3]*TM.wvx[n] + TM.rMatT[(n*9)+4]*TM.wvy[n] + TM.rMatT[(n*9)+5]*TM.wvz[n];
                 wvBz = TM.rMatT[(n*9)+6]*TM.wvx[n] + TM.rMatT[(n*9)+7]*TM.wvy[n] + TM.rMatT[(n*9)+8]*TM.wvz[n];
                 TM.rrvxf[n] = TM.mInertX[0]*wvBx*wvBx/TM.vvrSum;
                 TM.rrvyf[n] = TM.mInertY[0]*wvBy*wvBy/TM.vvrSum;
                 TM.rrvzf[n] = TM.mInertZ[0]*wvBx*wvBx/TM.vvrSum;
             }
         }else if(TM.ThermostatType == 5){
             TM.rotFrac = TM.rotationalTemperature/TM.instantTemperature;
             TM.tranFrac = TM.translationalTemperature/TM.instantTemperature;
             for(int n = 0; n < TM.nMol; n++){
                 TM.rvxf[n] = TM.rvx[n]/TM.vvSum;
                 TM.rvyf[n] = TM.rvy[n]/TM.vvSum;
                 TM.rvzf[n] = TM.rvz[n]/TM.vvSum;
                 TM.rrvxf[n] = TM.wvx[n]/TM.vvrSum;
                 TM.rrvyf[n] = TM.wvy[n]/TM.vvrSum;
                 TM.rrvzf[n] = TM.wvz[n]/TM.vvrSum;
             }
         }
    }
    public void AccumProps(TheMatrix TM, int icode){
        if(icode == 0){
            TM.totEnergySum = TM.kinEnergySum = TM.potEnergySum = 0.0;
            TM.potEnergyVdwSum = TM.potEnergyEeSum = TM.potEnergyRf1Sum = TM.potEnergyRf2Sum = TM.potEnergyEfSum = 0;
            TM.trzKinEnergySum = TM.rotKinEnergySum = 0.0;
            TM.aveTemperature = TM.pressure = TM.pressureSum = 0;
            TM.rotTempAve = TM.trzTempAve = 0;
        }else if(icode == 1){
            TM.totEnergySum = TM.totEnergySum + TM.totEnergyVal;
            TM.potEnergySum = TM.potEnergySum + TM.potEnergyVal;
            TM.potEnergyVdwSum = TM.potEnergyVdwSum + TM.potEnergyVdwVal;
            TM.potEnergyEeSum = TM.potEnergyEeSum + TM.potEnergyEeVal;
            TM.potEnergyRf1Sum = TM.potEnergyRf1Sum + TM.potEnergyRf1Val;
            TM.potEnergyRf2Sum = TM.potEnergyRf2Sum + TM.potEnergyRf2Val;
            TM.potEnergyEfSum = TM.potEnergyEfSum + TM.potEnergyEfVal;
            TM.kinEnergySum = TM.kinEnergySum + TM.kinEnergyVal;
            TM.trzKinEnergySum = TM.trzKinEnergySum + TM.trzKinEnergyVal;
            TM.rotKinEnergySum = TM.rotKinEnergySum + TM.rotKinEnergyVal;
            TM.aveTemperature = TM.aveTemperature + TM.instantTemperature;
            TM.rotTempAve = TM.rotTempAve + TM.rotationalTemperature;
            TM.trzTempAve = TM.trzTempAve + TM.translationalTemperature;
            TM.pressureSum = TM.pressureSum + TM.pressure;
        }else if(icode == 2){
            TM.totEnergySum /= TM.stepAvg;
            TM.potEnergySum /= TM.stepAvg;
            TM.potEnergyVdwSum /= TM.stepAvg;
            TM.potEnergyEeSum /= TM.stepAvg;
            TM.potEnergyRf1Sum /= TM.stepAvg;
            TM.potEnergyRf2Sum /= TM.stepAvg;
            TM.potEnergyEfSum /= TM.stepAvg;
            TM.kinEnergySum /= TM.stepAvg;
            TM.trzKinEnergySum /= TM.stepAvg;
            TM.rotKinEnergySum /= TM.stepAvg;
            TM.aveTemperature /= TM.stepAvg;
            TM.rotTempAve /= TM.stepAvg;
            TM.trzTempAve /= TM.stepAvg;
            TM.pressureSum /= TM.stepAvg;
        }
    }

    public void ComputeSiteForces(TheMatrix TM){
        double drx, dry, drz, shiftx, shifty, shiftz;
        double fcValEE, fcValRF, rr1, rr2, rr3, rrCut2, rrCut3, rri, rri3, uVal, uValVDW, uValEE;
        double fcValx, fcValy, fcValz, virSumTemp;
        int m1, m2, j1, j2, ms1F, ms1T, ms2F, ms2T, typeSum;
        rr2 = rri = 0;
        m1 = m2 = j1 = j2 = ms1F = ms1T = ms2F = ms2T = typeSum = 0;
        uValVDW = uValEE = fcValEE = fcValRF = 0;
        fcValx = fcValy = fcValz = virSumTemp = 0;
        rrCut2 = TM.rCut*TM.rCut;
        rrCut3 = rrCut2*TM.rCut;
        for(int n = 0; n < TM.totalSitesNumber; n++){
            TM.fxs[n] = 0.0; TM.fys[n] = 0.0;  TM.fzs[n] = 0.0;
        }
        TM.uSum = TM.uSumVDW = TM.uSumEE = TM.uSumRF1 = TM.uSumRF2 = TM.virSum = 0;
        for(m1 = 0; m1 < TM.nMol-1; m1++){
            ms1T = ms1T + TM.sitesMolIdx[m1];
            ms2F = ms1T;
            for(m2 = m1+1; m2 < TM.nMol; m2++){
                ms2T = ms2F + TM.sitesMolIdx[m2];
                drx = TM.rx[m1] - TM.rx[m2];           
                dry = TM.ry[m1] - TM.ry[m2];          
                drz = TM.rz[m1] - TM.rz[m2];
		shiftx = shifty = shiftz = 0.0;
		if(drx >= 0.5*TM.regionX)  shiftx = shiftx - TM.regionX;   
                else if(drx < -0.5*TM.regionX) shiftx = shiftx + TM.regionX; 
                if(dry >= 0.5*TM.regionY)  shifty = shifty - TM.regionY;
                else if(dry < -0.5*TM.regionY) shifty = shifty + TM.regionY;
                if(drz >= 0.5*TM.regionZ)  shiftz = shiftz - TM.regionZ;
                else if(drz < -0.5*TM.regionZ) shiftz = shiftz + TM.regionZ;
                drx = drx + shiftx;  dry = dry + shifty;  drz = drz + shiftz;
                rr2 = drx*drx + dry*dry + drz*drz;
                if(rr2 < rrCut2){
                    for(j1 = ms1F; j1 < ms1T; j1++){
                        for(j2 = ms2F; j2 < ms2T; j2++){
                            uValVDW = uValEE = fcValEE = fcValRF = 0;
                            typeSum =  TM.typeF[j1] + TM.typeF[j2];
                            if(TM.typeF[j1] == TM.typeF[j2] || typeSum == 5 || typeSum == 10){
                                drx = TM.rxs[j1] - TM.rxs[j2];
                                dry = TM.rys[j1] - TM.rys[j2];   
                                drz = TM.rzs[j1] - TM.rzs[j2];
                                drx = drx + shiftx;  dry = dry + shifty;  drz = drz + shiftz;
                                rr2 = drx*drx + dry*dry + drz*drz;
                                rr1 = Math.sqrt(rr2);
                                rr3 = rr2*rr1;                               
                                rri = 1.0/rr2;
                                // For now bcon is the one for the water model. This is becuase the LJ particles will use the same
                                // constant for now and until this procedure is adapted to deal with interaction between disimilar molecules.
                                switch(typeSum){
                                    case 2:
                                        rri3 = rri*rri*rri;
                                        uValVDW = 4 * rri3 * (rri3 - 1);
                                        fcValEE = 48 * rri3 * (rri3 - 0.5) * rri;
                                        break;
                                    case 4:
                                        uValEE =  4 * TM.bCon[0]/rr1;
                                        fcValEE = uValEE * rri;
                                        break;
                                    case 5:
                                        uValEE = -2 * TM.bCon[0]/rr1;
                                        fcValEE = uValEE * rri;
                                        break;
                                    case 6:                                       
                                        uValEE = TM.bCon[0]/rr1;
					fcValEE = uValEE * rri;
                                        break;
                                    case 12:
                                        // LJ particle is negative and repels oxygen
                                        uValEE = 2 * TM.LJcharge * TM.bCon[0]/rr1;
                                        fcValEE = uValEE * rri;
                                        break;
                                    case 13:
                                        // LJ particle is negative and attracts hydrogens.
                                        uValEE = -1 * TM.LJcharge * TM.bCon[0]/rr1;
                                        fcValEE = uValEE * rri;
                                        break;
                                    case 17:
                                        // LJ particle is negative and attracts hydrogens.
                                        uValEE = -1 * TM.LJcharge * TM.bCon[0]/rr1;
                                        fcValEE = uValEE * rri;
                                        break;
                                    case 10:
                                        uValEE = -1 * TM.bCon[0]/rr1;
                                        fcValEE = uValEE * rri;
                                        break;
                                    case 14:                    
                                        uValEE = TM.bCon[0]/rr1;
                                        fcValEE = uValEE * rri;
                                        break;
                                }
                                fcValx = (fcValEE)*drx;  fcValy = (fcValEE)*dry;     fcValz = (fcValEE)*drz;
                                TM.fxs[j1] = TM.fxs[j1] + fcValx;
                                TM.fys[j1] = TM.fys[j1] + fcValy;
                                      if(TM.typeF[j1] == 1){  TM.fzs[j1] = TM.fzs[j1] + fcValz;
                                }else if(TM.typeF[j1] == 2){  TM.fzs[j1] = TM.fzs[j1] + fcValz + (TM.electricField*(-2));
                                }else if(TM.typeF[j1] == 3){  TM.fzs[j1] = TM.fzs[j1] + fcValz + (TM.electricField);
                                }else if(TM.typeF[j1] == 7){  TM.fzs[j1] = TM.fzs[j1] + fcValz + (TM.electricField*(-1));}
                                TM.fxs[j2] = TM.fxs[j2] - fcValx;
                                TM.fys[j2] = TM.fys[j2] - fcValy;
                                      if(TM.typeF[j2] == 1){  TM.fzs[j2] = TM.fzs[j2] - fcValz;
                                }else if(TM.typeF[j2] == 2){  TM.fzs[j2] = TM.fzs[j2] - fcValz + (TM.electricField*(-2));
                                }else if(TM.typeF[j2] == 3){  TM.fzs[j2] = TM.fzs[j2] - fcValz + (TM.electricField);
                                }else if(TM.typeF[j2] == 7){  TM.fzs[j2] = TM.fzs[j2] - fcValz + (TM.electricField*(-1));}
                                TM.virSum = TM.virSum + (fcValx*drx + fcValy*dry + fcValz*drz); 
                                TM.uSumVDW += uValVDW;
                                TM.uSumEE += uValEE;
                                TM.uSum += (uValVDW + uValEE);
                                //System.out.println(m1+" "+m2+" "+j1+" "+j2+" vdw="+TM.uSumVDW+"+EE="+TM.uSumEE+"="+TM.uSum);
                            }
                        }
                    }
                 }
                 ms2F = ms2T;
            }
            ms1F = ms1T;
	}
    }
    public void ComputeSiteForcesRF(TheMatrix TM){
        double drx, dry, drz, shiftx, shifty, shiftz;
        double fcValEE, fcValRF, rr1, rr2, rr3, k1, k2, rrCut2, rrCut3, rri, rri3, uVal, uValVDW, uValEE, uValRF1;
        double fcValx, fcValy, fcValz, virSumTemp;
        int ms1F, ms1T, ms2F, ms2T, typeSum, Erf, Er;
        ms1F = ms1T = ms2F = ms2T = typeSum = 0;
        Erf = 78;    Er = 1;
        uValVDW = uValEE = fcValEE = fcValRF = uValRF1 = 0;
        fcValx = fcValy = fcValz = virSumTemp = 0;
        rrCut2 = TM.rCut*TM.rCut;
        rrCut3 = rrCut2*TM.rCut;
        k1 = (double) (Erf - Er)/(2*Erf+Er);  // (Erf - Er)/(2*Erf+Er)
        k2 = (double) (3*Erf)/(2*Erf+Er);     // (3*Erf)/(2*Erf+Er)
        for(int n = 0; n < TM.totalSitesNumber; n++){
              TM.fxs[n] = 0.0; TM.fys[n] = 0.0;  TM.fzs[n] = 0.0;
        }
        TM.uSum = TM.uSumVDW = TM.uSumEE = TM.uSumRF1 = TM.uSumRF2 = TM.virSum = 0;
        uValRF1 = -(2*TM.nMol*Math.PI*TM.density_convFact[0]*4)/(3*rrCut3);   // Cavity correction for when Erf = infinity.
        for(int m1 = 0; m1 < TM.nMol-1; m1++){
            ms1T = ms1T + TM.sitesMolIdx[m1];
            ms2F = ms1T;
            for(int m2 = m1+1; m2 < TM.nMol; m2++){
                ms2T = ms2F + TM.sitesMolIdx[m2];
                drx = TM.rx[m1] - TM.rx[m2];           
                dry = TM.ry[m1] - TM.ry[m2];          
                drz = TM.rz[m1] - TM.rz[m2];
		shiftx = shifty = shiftz = 0.0;
		if(drx >= 0.5*TM.regionX)  shiftx = shiftx - TM.regionX; 
                else if(drx < -0.5*TM.regionX) shiftx = shiftx + TM.regionX; 
                if(dry >= 0.5*TM.regionY)  shifty = shifty - TM.regionY;
                else if(dry < -0.5*TM.regionY) shifty = shifty + TM.regionY;
                if(drz >= 0.5*TM.regionZ)  shiftz = shiftz - TM.regionZ;
                else if(drz < -0.5*TM.regionZ) shiftz = shiftz + TM.regionZ;
                drx = drx + shiftx;  dry = dry + shifty;  drz = drz + shiftz;
                rr2 = drx*drx + dry*dry + drz*drz;  
                if(rr2 < rrCut2){
                    for(int j1 = ms1F; j1 < ms1T; j1++){
                        for(int j2 = ms2F; j2 < ms2T; j2++){
                            uValVDW = uValEE = fcValEE = fcValRF = 0;
                            typeSum =  TM.typeF[j1] + TM.typeF[j2];
                            if(TM.typeF[j1] == TM.typeF[j2] || typeSum == 5 || typeSum == 10){
                                drx = TM.rxs[j1] - TM.rxs[j2];
                                dry = TM.rys[j1] - TM.rys[j2];   
                                drz = TM.rzs[j1] - TM.rzs[j2];
                                drx = drx + shiftx;  dry = dry + shifty;  drz = drz + shiftz;
                                rr2 = drx*drx + dry*dry + drz*drz;
                                rr1 = Math.sqrt(rr2); 
                                rr3 = rr2*rr1;                               
                                rri = 1.0/rr2;
                                switch(typeSum){
                                    case 2:
                                        rri3 = rri*rri*rri;
                                        uValVDW = 4 * rri3 * (rri3 - 1);
                                        fcValEE = 48 * rri3 * (rri3 - 0.5) * rri;
                                        break;
                                    case 4:
                                        uValEE =  4 * TM.bCon[0]*((1/rr1) + (k1 * rr2/rrCut3) - (k2/TM.rCut));
                                        fcValEE = uValEE * rri;
                                        break;
                                    case 5:
                                        uValEE = -2 * TM.bCon[0]*((1/rr1) + (k1 * rr2/rrCut3) - (k2/TM.rCut));
                                        fcValEE = uValEE * rri;
                                        break;
                                    case 6:                                       
                                        uValEE = TM.bCon[0]*((1/rr1) + (k1 * rr2/rrCut3) - (k2/TM.rCut));
					fcValEE = uValEE * rri;
                                        break;
                                    case 12:
                                        // LJ particle is negative and repels oxygen
                                        uValEE = 2 * TM.LJcharge * TM.bCon[0]/rr1;
                                        fcValEE = uValEE * rri;
                                        break;
                                    case 13:
                                        // LJ particle is negative and attracts hydrogens.
                                        uValEE = -1 * TM.LJcharge * TM.bCon[0]/rr1;
                                        fcValEE = uValEE * rri;
                                        break;
                                    case 17:
                                        // LJ particle is negative and attracts hydrogens.
                                        uValEE = -1 * TM.LJcharge * TM.bCon[0]/rr1;
                                        fcValEE = uValEE * rri;
                                        break;
                                    case 10:                                        
                                        uValEE = -1 * TM.bCon[0]*((1/rr1) + (k1 * rr2/rrCut3) - (k2/TM.rCut));
                                        fcValEE = uValEE * rri;
                                        break;
                                    case 14:                                      
                                        uValEE = TM.bCon[0]*((1/rr1) + (k1 * rr2/rrCut3) - (k2/TM.rCut));
                                        fcValEE = uValEE * rri;
                                        break;
                                }
                                fcValx = (fcValEE)*drx;  fcValy = (fcValEE)*dry;     fcValz = (fcValEE)*drz;

                                TM.fxs[j1] = TM.fxs[j1] + fcValx;
                                TM.fys[j1] = TM.fys[j1] + fcValy;
                                      if(TM.typeF[j1] == 1){  TM.fzs[j1] = TM.fzs[j1] + fcValz;
                                }else if(TM.typeF[j1] == 2){  TM.fzs[j1] = TM.fzs[j1] + fcValz + (TM.electricField*(-2));
                                }else if(TM.typeF[j1] == 3){  TM.fzs[j1] = TM.fzs[j1] + fcValz + (TM.electricField);
                                }else if(TM.typeF[j1] == 7){  TM.fzs[j1] = TM.fzs[j1] + fcValz + (TM.electricField*(-1));}
                                
                                TM.fxs[j2] = TM.fxs[j2] - fcValx;
                                TM.fys[j2] = TM.fys[j2] - fcValy;
                                      if(TM.typeF[j2] == 1){  TM.fzs[j2] = TM.fzs[j2] - fcValz;
                                }else if(TM.typeF[j2] == 2){  TM.fzs[j2] = TM.fzs[j2] - fcValz + (TM.electricField*(-2));
                                }else if(TM.typeF[j2] == 3){  TM.fzs[j2] = TM.fzs[j2] - fcValz + (TM.electricField);
                                }else if(TM.typeF[j2] == 7){  TM.fzs[j2] = TM.fzs[j2] - fcValz + (TM.electricField*(-1));}
                                TM.virSum = TM.virSum + (fcValx*drx + fcValy*dry + fcValz*drz); 
                                TM.uSumVDW += uValVDW;
                                TM.uSumEE += uValEE;
                                TM.uSum += (uValVDW + uValEE);
                            }
                        }
                    }
                 }
                 ms2F = ms2T;
            }
            ms1F = ms1T;
	}
        TM.uSumRF1 = uValRF1;
	TM.uSum += TM.uSumRF1;
    }
    public void ComputeSiteForcesWxyz(TheMatrix TM){
        double drx, dry, drz;
        double fcValEE, fcValRF, rr1, rr2, rr3, rrCut2, rrCut3, rri, rri3, uVal, uValVDW, uValEE;
        double fcValx, fcValy, fcValz, virSumTemp;
        int ms1F, ms1T, ms2F, ms2T, typeSum;
        ms1F = ms1T = ms2F = ms2T = typeSum = 0;
        uValVDW = uValEE = fcValEE = fcValRF = 0;
        fcValx = fcValy = fcValz = virSumTemp = 0;
        rrCut2 = TM.rCut*TM.rCut;
        rrCut3 = rrCut2*TM.rCut;
        for(int n = 0; n < TM.totalSitesNumber; n++){
              TM.fxs[n] = 0.0; TM.fys[n] = 0.0;  TM.fzs[n] = 0.0;
        }
        TM.uSum = TM.uSumVDW = TM.uSumEE = TM.uSumRF1 = TM.uSumRF2 = TM.virSum = 0;
        for(int m1 = 0; m1 < TM.nMol-1; m1++){
            ms1T = ms1T + TM.sitesMolIdx[m1];
            ms2F = ms1T;
            for(int m2 = m1+1; m2 < TM.nMol; m2++){
                ms2T = ms2F + TM.sitesMolIdx[m2];
                drx = TM.rx[m1] - TM.rx[m2];           
                dry = TM.ry[m1] - TM.ry[m2];          
                drz = TM.rz[m1] - TM.rz[m2];
                rr2 = drx*drx + dry*dry + drz*drz;  
                if(rr2 < rrCut2){
                    for(int j1 = ms1F; j1 < ms1T; j1++){
                        for(int j2 = ms2F; j2 < ms2T; j2++){
                            uValVDW = uValEE = fcValEE = fcValRF = 0;
                            typeSum =  TM.typeF[j1] + TM.typeF[j2];
                            if(TM.typeF[j1] == TM.typeF[j2] || typeSum == 5 || typeSum == 10){
                                drx = TM.rxs[j1] - TM.rxs[j2];
                                dry = TM.rys[j1] - TM.rys[j2];   
                                drz = TM.rzs[j1] - TM.rzs[j2];
                                rr2 = drx*drx + dry*dry + drz*drz;
                                rr1 = Math.sqrt(rr2); 
                                rr3 = rr2*rr1;                               
                                rri = 1.0/rr2;
                                switch(typeSum){
                                    case 2:
                                        rri3 = rri*rri*rri;
                                        uValVDW = 4 * rri3 * (rri3 - 1);
                                        fcValEE = 48 * rri3 * (rri3 - 0.5) * rri;
                                        break;
                                    case 4:
                                        uValEE =  4 * TM.bCon[0]/rr1;
                                        fcValEE = uValEE * rri;
                                        break;
                                    case 5:
                                        uValEE = -2 * TM.bCon[0]/rr1;
                                        fcValEE = uValEE * rri;
                                        break;
                                    case 6:                                       
                                        uValEE = TM.bCon[0]/rr1;
					fcValEE = uValEE * rri;
                                        break;
                                    case 12:
                                        // LJ particle is negative and repels oxygen
                                        uValEE = 2 * TM.LJcharge * TM.bCon[0]/rr1;
                                        fcValEE = uValEE * rri;
                                        break;
                                    case 13:
                                        // LJ particle is negative and attracts hydrogens.
                                        uValEE = -1 * TM.LJcharge * TM.bCon[0]/rr1;
                                        fcValEE = uValEE * rri;
                                        break;
                                    case 17:
                                        // LJ particle is negative and attracts hydrogens.
                                        uValEE = -1 * TM.LJcharge * TM.bCon[0]/rr1;
                                        fcValEE = uValEE * rri;
                                        break;
                                    case 10:                                        
                                        uValEE = -1 * TM.bCon[0]/rr1;
                                        fcValEE = uValEE * rri;
                                        break;
                                    case 14:                                        
                                        uValEE = TM.bCon[0]/rr1;
                                        fcValEE = uValEE * rri;
                                        break;
                                }
                                fcValx = (fcValEE)*drx;  fcValy = (fcValEE)*dry;     fcValz = (fcValEE)*drz;

                                TM.fxs[j1] = TM.fxs[j1] + fcValx;
                                TM.fys[j1] = TM.fys[j1] + fcValy;
                                      if(TM.typeF[j1] == 1){  TM.fzs[j1] = TM.fzs[j1] + fcValz;
                                }else if(TM.typeF[j1] == 2){  TM.fzs[j1] = TM.fzs[j1] + fcValz + (TM.electricField*(-2));
                                }else if(TM.typeF[j1] == 3){  TM.fzs[j1] = TM.fzs[j1] + fcValz + (TM.electricField);
                                }else if(TM.typeF[j1] == 7){  TM.fzs[j1] = TM.fzs[j1] + fcValz + (TM.electricField*(-1));}
                                
                                TM.fxs[j2] = TM.fxs[j2] - fcValx;
                                TM.fys[j2] = TM.fys[j2] - fcValy;
                                      if(TM.typeF[j2] == 1){  TM.fzs[j2] = TM.fzs[j2] - fcValz;
                                }else if(TM.typeF[j2] == 2){  TM.fzs[j2] = TM.fzs[j2] - fcValz + (TM.electricField*(-2));
                                }else if(TM.typeF[j2] == 3){  TM.fzs[j2] = TM.fzs[j2] - fcValz + (TM.electricField);
                                }else if(TM.typeF[j2] == 7){  TM.fzs[j2] = TM.fzs[j2] - fcValz + (TM.electricField*(-1));}
                                TM.virSum = TM.virSum + (fcValx*drx + fcValy*dry + fcValz*drz); 
                                TM.uSumVDW += uValVDW;
                                TM.uSumEE += uValEE;
                                TM.uSum += (uValVDW + uValEE);
                            }
                        }
                    }
                }
                ms2F = ms2T;
            }
            ms1F = ms1T;
	}
    }
    public void ComputeSiteForcesWz(TheMatrix TM){
        double drx, dry, drz, shiftx, shifty;// shiftz;
        double fcValEE, fcValRF, rr1, rr2, rr3, rrCut2, rrCut3, rri, rri3, uVal, uValVDW, uValEE;
        double fcValx, fcValy, fcValz, virSumTemp;
        int ms1F, ms1T, ms2F, ms2T, typeSum;
        ms1F = ms1T = ms2F = ms2T = typeSum = 0;
        uValVDW = uValEE = fcValEE = fcValRF = 0;
        fcValx = fcValy = fcValz = virSumTemp = 0;
        rrCut2 = TM.rCut*TM.rCut;
        rrCut3 = rrCut2*TM.rCut;
        for(int n = 0; n < TM.totalSitesNumber; n++){
              TM.fxs[n] = 0.0; TM.fys[n] = 0.0;  TM.fzs[n] = 0.0;
        }
        TM.uSum = TM.uSumVDW = TM.uSumEE = TM.uSumRF1 = TM.uSumRF2 = TM.virSum = 0;
        for(int m1 = 0; m1 < TM.nMol-1; m1++){
            ms1T = ms1T + TM.sitesMolIdx[m1];
            ms2F = ms1T;
            for(int m2 = m1+1; m2 < TM.nMol; m2++){
                ms2T = ms2F + TM.sitesMolIdx[m2];
                drx = TM.rx[m1] - TM.rx[m2];           
                dry = TM.ry[m1] - TM.ry[m2];          
                drz = TM.rz[m1] - TM.rz[m2];
		shiftx = shifty = 0.0; //shiftz = 0.0;
		if(drx >= 0.5*TM.regionX)  shiftx = shiftx - TM.regionX;  
                else if(drx < -0.5*TM.regionX) shiftx = shiftx + TM.regionX; 
                if(dry >= 0.5*TM.regionY)  shifty = shifty - TM.regionY;
                else if(dry < -0.5*TM.regionY) shifty = shifty + TM.regionY;
        //      if(drz >= 0.5*TM.regionZ)  shiftz = shiftz - TM.regionZ;
        //      else if(drz < -0.5*TM.regionZ) shiftz = shiftz + TM.regionZ;
                drx = drx + shiftx;  dry = dry + shifty;  //drz = drz + shiftz;
                rr2 = drx*drx + dry*dry + drz*drz;  
                if(rr2 < rrCut2){
                    for(int j1 = ms1F; j1 < ms1T; j1++){
                        for(int j2 = ms2F; j2 < ms2T; j2++){
                            uValVDW = uValEE = fcValEE = fcValRF = 0;
                            typeSum =  TM.typeF[j1] + TM.typeF[j2];
                            if(TM.typeF[j1] == TM.typeF[j2] || typeSum == 5 || typeSum == 10){
                                drx = TM.rxs[j1] - TM.rxs[j2];
                                dry = TM.rys[j1] - TM.rys[j2];   
                                drz = TM.rzs[j1] - TM.rzs[j2];
                                drx = drx + shiftx;  dry = dry + shifty;  //drz = drz + shiftz;
                                rr2 = drx*drx + dry*dry + drz*drz;
                                rr1 = Math.sqrt(rr2); 
                                rr3 = rr2*rr1;                               
                                rri = 1.0/rr2;
                                switch(typeSum){
                                    case 2:
                                        rri3 = rri*rri*rri;
                                        uValVDW = 4 * rri3 * (rri3 - 1);
                                        fcValEE = 48 * rri3 * (rri3 - 0.5) * rri;
                                        break;
                                    case 4:
                                        uValEE =  4 * TM.bCon[0]/rr1;
                                        fcValEE = uValEE * rri;
                                        break;
                                    case 5:
                                        uValEE = -2 * TM.bCon[0]/rr1;
                                        fcValEE = uValEE * rri;
                                        break;
                                    case 6:                                       
                                        uValEE = TM.bCon[0]/rr1;
					fcValEE = uValEE * rri;
                                        break;
                                    case 12:
                                        // LJ particle is negative and repels oxygen
                                        uValEE = 2 * TM.LJcharge * TM.bCon[0]/rr1;
                                        fcValEE = uValEE * rri;
                                        break;
                                    case 13:
                                        // LJ particle is negative and attracts hydrogens.
                                        uValEE = -1 * TM.LJcharge * TM.bCon[0]/rr1;
                                        fcValEE = uValEE * rri;
                                        break;
                                    case 17:
                                        // LJ particle is negative and attracts hydrogens.
                                        uValEE = -1 * TM.LJcharge * TM.bCon[0]/rr1;
                                        fcValEE = uValEE * rri;
                                        break;
                                    case 10:                                        
                                        uValEE = -1 * TM.bCon[0]/rr1;
                                        fcValEE = uValEE * rri;
                                        break;
                                    case 14:                                        
                                        uValEE = TM.bCon[0]/rr1;
                                        fcValEE = uValEE * rri;
                                        break;
                                }
                                fcValx = (fcValEE)*drx;  fcValy = (fcValEE)*dry;     fcValz = (fcValEE)*drz;

                                TM.fxs[j1] = TM.fxs[j1] + fcValx;
                                TM.fys[j1] = TM.fys[j1] + fcValy;
                                      if(TM.typeF[j1] == 1){  TM.fzs[j1] = TM.fzs[j1] + fcValz;
                                }else if(TM.typeF[j1] == 2){  TM.fzs[j1] = TM.fzs[j1] + fcValz + (TM.electricField*(-2));
                                }else if(TM.typeF[j1] == 3){  TM.fzs[j1] = TM.fzs[j1] + fcValz + (TM.electricField);
                                }else if(TM.typeF[j1] == 7){  TM.fzs[j1] = TM.fzs[j1] + fcValz + (TM.electricField*(-1));}
                                
                                TM.fxs[j2] = TM.fxs[j2] - fcValx;
                                TM.fys[j2] = TM.fys[j2] - fcValy;
                                      if(TM.typeF[j2] == 1){ TM.fzs[j2] = TM.fzs[j2] - fcValz;
                                }else if(TM.typeF[j2] == 2){ TM.fzs[j2] = TM.fzs[j2] - fcValz + (TM.electricField*(-2));
                                }else if(TM.typeF[j2] == 3){ TM.fzs[j2] = TM.fzs[j2] - fcValz + (TM.electricField);
				}else if(TM.typeF[j2] == 7){ TM.fzs[j2] = TM.fzs[j2] - fcValz + (TM.electricField*(-1));}
                                TM.virSum = TM.virSum + (fcValx*drx + fcValy*dry + fcValz*drz); 
                                TM.uSumVDW += uValVDW;
                                TM.uSumEE += uValEE;
                                TM.uSum += (uValVDW + uValEE);
                            }
                        }
                    }
                }
                ms2F = ms2T;
            }
            ms1F = ms1T;
	}
    }
}
