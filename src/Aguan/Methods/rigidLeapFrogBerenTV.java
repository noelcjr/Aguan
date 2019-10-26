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
public class rigidLeapFrogBerenTV extends MD{
    public rigidLeapFrogBerenTV(String[] args){
        super(args);
        AccumProps(TM, 0);
             if(TM.boundaryCondition.equals("periodic")){    super.ComputeSiteForces(TM);}
        else if(TM.boundaryCondition.equals("xyzWall")){      ComputeSiteForcesWxyz(TM);  ComputeWallForcesXYZ(TM);}
        else if(TM.boundaryCondition.equals("zWall")){        ComputeSiteForcesWz(TM);    ComputeWallForcesZ(TM);}
        ComputeTorqs(TM);
        EvalProps(TM);
        AccumProps(TM, 1);
        System.out.println("stepCount = "+TM.stepCount);
        if(TM.stepCount == 0){
           timeNow = TM.stepCount * TM.deltaT;
           AdjustTemp(TM);
           OUT.PrintSummary(TM, timeNow);
           RR.writeRestart(TM);
           AccumProps(TM,0);
        }
        if(TM.boundaryCondition.equals("periodic")){
            while(TM.moreCycles){
                SingleStep(TM,RR);
                if(TM.stepCount >= TM.stepLimit){
                   TM.moreCycles = false;
                }
            }
        }else if(TM.boundaryCondition.equals("xyzWall")){
            while(TM.moreCycles){
                SingleStepXYZ(TM,RR);
                if(TM.stepCount >= TM.stepLimit){
                   TM.moreCycles = false;
                }
            }
        }else if(TM.boundaryCondition.equals("zWall")){
            while(TM.moreCycles){
                SingleStepZ(TM,RR);
                if(TM.stepCount >= TM.stepLimit){
                   TM.moreCycles = false;
                }
            }
        }else{
             System.out.println("INPUT ERROR: value for boundary condition is not valid, only");
             System.out.println("the case sensitive values of periodic, xyzWall and zWall are allowed.");
        }
        try{
        dcdM.data_out.close();
        }catch( IOException e ){System.err.println( e );}
    }
    public void SingleStepXYZ(TheMatrix TM, restart RR){
        TM.stepCount++;
        timeNow = TM.stepCount * TM.deltaT;
        TM.targetTemperature = TM.targetTemperature + TM.TempRateChange;    
        LeapfrogStep(1, TM);
        GenSiteCoord(TM); 
        ComputeSiteForcesWxyz(TM);
        ComputeWallForcesXYZ(TM);
        ComputeTorqs(TM);
        AdjustTemp(TM);
        LeapfrogStep(2, TM);
        ApplyBoundaryCond(TM);
        EvalProps(TM);
        if(TM.stepAvg > 1){ 
            AccumProps(TM, 1); 
            if(TM.stepCount % TM.stepAvg == 0){ 
               AccumProps(TM, 2); 
               OUT.PrintSummary_ave(TM,timeNow);
               AccumProps(TM, 0); 
            }   
        }else{
            AccumProps(TM, 1); 
            OUT.PrintSummary(TM,timeNow);
            AccumProps(TM, 0);
        }
        if(TM.stepCount % TM.trajOut == 0){
            dcdM.write_dcdStep(TM.sitesNoVDW,TM.nMol,TM.atomType,TM.sitesMolIdx,TM.rxs,TM.rys,TM.rzs,TM.ro[0]);
            RR.writeRestart(TM);
        }
    }
    public void SingleStepZ(TheMatrix TM, restart RR){
        TM.stepCount++;
        timeNow = TM.stepCount * TM.deltaT;
        TM.targetTemperature = TM.targetTemperature + TM.TempRateChange;    
        LeapfrogStep(1, TM);
        GenSiteCoord(TM); 
        ComputeSiteForcesWz(TM);
        ComputeWallForcesZ(TM);
        ComputeTorqs(TM);
        AdjustTemp(TM);
        LeapfrogStep(2, TM);
        ApplyBoundaryCond(TM);
        EvalProps(TM);
        if(TM.stepAvg > 1){ 
            AccumProps(TM, 1); 
            if(TM.stepCount % TM.stepAvg == 0){ 
               AccumProps(TM, 2); 
               OUT.PrintSummary_ave(TM,timeNow);
               AccumProps(TM, 0); 
            }   
        }else{
            AccumProps(TM, 1); 
            OUT.PrintSummary(TM,timeNow);
            AccumProps(TM, 0);
        }
        if(TM.stepCount % TM.trajOut == 0){
            dcdM.write_dcdStep(TM.sitesNoVDW,TM.nMol,TM.atomType,TM.sitesMolIdx,TM.rxs,TM.rys,TM.rzs,TM.ro[0]);
            RR.writeRestart(TM);
        }
    }
    public void SingleStep(TheMatrix TM, restart RR){
        TM.stepCount++;
        timeNow = TM.stepCount * TM.deltaT;
        TM.targetTemperature = TM.targetTemperature + TM.TempRateChange;    
        LeapfrogStep(1, TM);
        GenSiteCoord(TM); 
        ComputeSiteForces(TM);
        ComputeTorqs(TM);
        AdjustTemp(TM);
        LeapfrogStep(2, TM);
        ApplyBoundaryCond(TM);
        EvalProps(TM);
        if(TM.stepAvg > 1){
            AccumProps(TM, 1);
            if(TM.stepCount % TM.stepAvg == 0){
               AccumProps(TM, 2);
               OUT.PrintSummary_ave(TM,timeNow);
               AccumProps(TM, 0);
            }
        }else{
            AccumProps(TM, 1);
            OUT.PrintSummary(TM,timeNow);
            AccumProps(TM, 0);
        }
        if(TM.stepCount % TM.trajOut == 0){
            dcdM.write_dcdStep(TM.sitesNoVDW,TM.nMol,TM.atomType,TM.sitesMolIdx,TM.rxs,TM.rys,TM.rzs,TM.ro[0]);
            RR.writeRestart(TM);
        }
    }
    public void LeapfrogStep(int part, TheMatrix TM){
          double[] mc, mt;
          double tx, ty, tz;
          mc = new double[9];          mt = new double[9];
          double temp = 0.5*TM.deltaT;
          if(part == 1){
             for(int n = 0; n < TM.nMol; n++){
                 TM.wvx[n] = TM.wvx[n] + temp * TM.wax[n];
                 TM.wvy[n] = TM.wvy[n] + temp * TM.way[n];
                 TM.wvz[n] = TM.wvz[n] + temp * TM.waz[n];
                 TM.rvx[n] = TM.rvx[n] + temp * TM.rax[n];
                 TM.rvy[n] = TM.rvy[n] + temp * TM.ray[n];
                 TM.rvz[n] = TM.rvz[n] + temp * TM.raz[n];
             }
             for(int n = 0; n < TM.nMol; n++){
                 tx = temp * TM.wvx[n];
                 ty = temp * TM.wvy[n];
                 tz = temp * TM.wvz[n];
                 BuildStepRmat(mc, tx, ty, tz);
                 MulMat(mt, mc, TM.rMatT, 3, n);
                 TM.rMatT[n*9] = mt[0];      TM.rMatT[n*9+1] = mt[1];    TM.rMatT[n*9+2] = mt[2];
                 TM.rMatT[n*9+3] = mt[3];    TM.rMatT[n*9+4] = mt[4];    TM.rMatT[n*9+5] = mt[5];
                 TM.rMatT[n*9+6] = mt[6];    TM.rMatT[n*9+7] = mt[7];    TM.rMatT[n*9+8] = mt[8];
             }
             for(int n = 0; n < TM.nMol; n++){
                 TM.rx[n] = TM.rx[n] + TM.deltaT * TM.rvx[n];
                 TM.ry[n] = TM.ry[n] + TM.deltaT * TM.rvy[n];
                 TM.rz[n] = TM.rz[n] + TM.deltaT * TM.rvz[n];
             }
          }else{
             for(int n = 0; n < TM.nMol; n++){
                 TM.wvx[n] = TM.wvx[n] + temp * TM.wax[n];
                 TM.wvy[n] = TM.wvy[n] + temp * TM.way[n];
                 TM.wvz[n] = TM.wvz[n] + temp * TM.waz[n];
                 TM.rvx[n] = TM.rvx[n] + temp * TM.rax[n];
                 TM.rvy[n] = TM.rvy[n] + temp * TM.ray[n];
                 TM.rvz[n] = TM.rvz[n] + temp * TM.raz[n];
             }
          }
    }
    public void BuildStepRmat(double[] mp, double ax, double ay, double az){
        double[] m1 = new double[9];        double[] m2 = new double[9];
        double[] c = new double[3];         double[] s = new double[3];
        double ak, c0c2, c0s2, s0c2, s0s2, t;
        ak = 0;
        //#define VComp(v, k)   *((k == 0) ? &(v).x : ((k == 1) ? &(v).y : &(v).z))
        for(int k = 0; k < 3; k++){
            if(k == 0){  ak = ax;
            }else if(k == 1){   ak = ay;
            }else if(k == 2){   ak = az;}
            t = 0.25 * ak*ak;
            c[k] = (1.0-t)/(1.0+t);
            s[k] = ak/(1.0+t);
        }
        c0c2 = c[0] * c[2];
        c0s2 = c[0] * s[2];
        s0c2 = s[0] * c[2];
        s0s2 = s[0] * s[2];
        m1[0] = c[1] * c[2];
        m1[1] = s0c2 * s[1] + c0s2;
        m1[2] = -c0c2 * s[1] + s0s2;
        m1[3] = -c[1] * s[2];
        m1[4] = -s0s2 * s[1] + c0c2;
        m1[5] = c0s2 * s[1] + s0c2;
        m1[6] = s[1];
        m1[7] = -s[0] * c[1];
        m1[8] = c[0] * c[1];
        m2[0] = m1[0];
        m2[1] = -m1[3];
        m2[2] = -m1[6];
        m2[3] = s0c2 * s[1] - c0s2;
        m2[4] = s0s2 * s[1] + c0c2;
        m2[5] = -m1[7];
        m2[6] = c0c2 * s[1] + s0s2;
        m2[7] = c0s2 * s[1] - s0c2;
        m2[8] = m1[8];
        MulMat(mp, m1, m2, 3, 0);
    }
    public void MulMat(double[] a, double[] b, double[] c, int nn, int n){
      //  #define MAT(a, n, i, j)  (a)[(i) + n * (j)]
        for(int i = 0; i < nn; i++){
            for(int j = 0; j < nn; j++){
                a[i + nn * j] = 0;
                for(int k = 0; k < nn; k++){
                    a[i + (nn * j)] += b[i + (nn * k)] * c[(n*9)+(k + (nn * j))];
                }
            }
        }
    }
    public void ComputeTorqs(TheMatrix TM){
        double drx, dry, drz, tx, ty, tz, torqSx, torqSy, torqSz, waBx, waBy, waBz;
        int index = 0;
        for(int n = 0; n < TM.nMol; n++){
            TM.rax[n] = TM.ray[n] = TM.raz[n] = 0.0;
            torqSx = torqSy = torqSz = 0.0;
            for(int j = 0; j < TM.sitesMolIdx[n]; j++){
                TM.rax[n] = TM.rax[n] + TM.fxs[index + j];                   // Right assumption?  Translational accelerations are summed,
                TM.ray[n] = TM.ray[n] + TM.fys[index + j];                   // Rotational accelerations are then obtained from this
                TM.raz[n] = TM.raz[n] + TM.fzs[index + j];                   // translational accelerations. Shouldn't accelerations
                drx = TM.rxs[index + j] - TM.rx[n];                          // be distributed between tranlation and rotation? Instead,
                dry = TM.rys[index + j] - TM.ry[n];                          // what seems to occur is that all the accelerations are added
                drz = TM.rzs[index + j] - TM.rz[n];                          // to translation and then a fraction of those are used for rotational
                tx = dry*TM.fzs[index + j] - drz*TM.fys[index + j];          // accelerations. There seems to be a double count.
                ty = drz*TM.fxs[index + j] - drx*TM.fzs[index + j];          // some of the acceleration that goes to rotation
                tz = drx*TM.fys[index + j] - dry*TM.fxs[index + j];          // is already going to translation. Also, the rotational
                torqSx =  torqSx + tx;     torqSy =  torqSy + ty;     torqSz =  torqSz + tz; // inertia tensor is calculated based on masses only
            }                                                                                  // not charges. is this right?  could this be fixed?
            waBx = TM.rMatT[(n*9)+0]*torqSx + TM.rMatT[(n*9)+1]*torqSy + TM.rMatT[(n*9)+2]*torqSz;
            waBy = TM.rMatT[(n*9)+3]*torqSx + TM.rMatT[(n*9)+4]*torqSy + TM.rMatT[(n*9)+5]*torqSz;
            waBz = TM.rMatT[(n*9)+6]*torqSx + TM.rMatT[(n*9)+7]*torqSy + TM.rMatT[(n*9)+8]*torqSz;
            waBx /= TM.mInertX[0];
            waBy /= TM.mInertY[0];
            waBz /= TM.mInertZ[0];
            TM.wax[n] = TM.rMatT[(n*9)+0]*waBx + TM.rMatT[(n*9)+3]*waBy + TM.rMatT[(n*9)+6]*waBz;
            TM.way[n] = TM.rMatT[(n*9)+1]*waBx + TM.rMatT[(n*9)+4]*waBy + TM.rMatT[(n*9)+7]*waBz;
            TM.waz[n] = TM.rMatT[(n*9)+2]*waBx + TM.rMatT[(n*9)+5]*waBy + TM.rMatT[(n*9)+8]*waBz;
            index = index + TM.sitesMolIdx[n];
        }
    }
    public void AdjustTemp(TheMatrix TM){
        double wvBx, wvBy, wvBz;
        double vFac, vFacT, vFacR;
        //   vFacT = Math.sqrt(1+(TM.deltaT/TM.tau)*((TM.targetTemperature/(TM.instantTemperature*TM.tranFrac)) - 1));
        //   vFacR = Math.sqrt(1+(TM.deltaT/TM.tau)*((TM.targetTemperature/(TM.instantTemperature*TM.ro[0]tFrac)) - 1));
        if(TM.ThermostatType == 1){  // According with the Equipartition, No Berensen.
               // Temperture gradient cannot work with this thermostat unless velMag is updated
               // everytime the temperature is increased.
               vFacT = TM.velMag/Math.sqrt(TM.vvSum/(TM.nMol/2));
               TM.vvSum = 0.0;
               for(int n = 0; n < TM.nMol; n++){
                   TM.rvx[n] *= vFacT;     TM.rvy[n] *= vFacT;    TM.rvz[n] *= vFacT;
                   TM.vvSum += TM.rvx[n]*TM.rvx[n] + TM.rvy[n]*TM.rvy[n] + TM.rvz[n]*TM.rvz[n];
               }
               TM.translationalTemperature = 2*(TM.vvSum/TM.nMol)/3.0;
               TM.trzKinEnergyVal = 0.5 * TM.vvSum/TM.nMol;
               vFacR = TM.velMag/Math.sqrt(TM.vvrSum/(TM.nMol/2));
               TM.vvrSum = 0.0;
               for(int n = 0; n < TM.nMol; n++){
                   TM.wvx[n] *= vFacR;     TM.wvy[n] *= vFacR;    TM.wvz[n] *= vFacR;
                   wvBx = TM.rMatT[(n*9)+0]*TM.wvx[n] + TM.rMatT[(n*9)+1]*TM.wvy[n] + TM.rMatT[(n*9)+2]*TM.wvz[n];
                   wvBy = TM.rMatT[(n*9)+3]*TM.wvx[n] + TM.rMatT[(n*9)+4]*TM.wvy[n] + TM.rMatT[(n*9)+5]*TM.wvz[n];
                   wvBz = TM.rMatT[(n*9)+6]*TM.wvx[n] + TM.rMatT[(n*9)+7]*TM.wvy[n] + TM.rMatT[(n*9)+8]*TM.wvz[n];
                   TM.vvrSum += TM.mInertX[0]*wvBx*wvBx + TM.mInertY[0]*wvBy*wvBy + TM.mInertZ[0]*wvBz*wvBz;
               }
               TM.rotationalTemperature = 2*(TM.vvrSum/TM.nMol)/3.0;
               TM.rotKinEnergyVal = 0.5 * TM.vvrSum/TM.nMol;
               //System.out.println("TM.velMag="+TM.velMag+" vFacT="+vFacT+" vFacR="+vFacR+" vvSum="+TM.vvSum+" vvrSum="+TM.vvrSum);
        }else if(TM.ThermostatType == 2){
               vFacT = Math.sqrt(1+(TM.tau*(((0.5*TM.targetTemperature)/TM.translationalTemperature)-1)));
               TM.vvSum = 0.0;
               for(int n = 0; n < TM.nMol; n++){
                   TM.rvx[n] *= vFacT;     TM.rvy[n] *= vFacT;    TM.rvz[n] *= vFacT;
                   TM.vvSum += TM.rvx[n]*TM.rvx[n] + TM.rvy[n]*TM.rvy[n] + TM.rvz[n]*TM.rvz[n];
               }
               TM.translationalTemperature = (TM.vvSum/TM.nMol)/3.0;
               TM.trzKinEnergyVal = 0.5 * TM.vvSum/TM.nMol;
               vFacR = Math.sqrt(1+(TM.tau*(((0.5*TM.targetTemperature)/TM.rotationalTemperature)-1)));
               TM.vvrSum = 0.0;
               for(int n = 0; n < TM.nMol; n++){
                   TM.wvx[n] *= vFacR;     TM.wvy[n] *= vFacR;    TM.wvz[n] *= vFacR;
                   wvBx = TM.rMatT[(n*9)+0]*TM.wvx[n] + TM.rMatT[(n*9)+1]*TM.wvy[n] + TM.rMatT[(n*9)+2]*TM.wvz[n];
                   wvBy = TM.rMatT[(n*9)+3]*TM.wvx[n] + TM.rMatT[(n*9)+4]*TM.wvy[n] + TM.rMatT[(n*9)+5]*TM.wvz[n];
                   wvBz = TM.rMatT[(n*9)+6]*TM.wvx[n] + TM.rMatT[(n*9)+7]*TM.wvy[n] + TM.rMatT[(n*9)+8]*TM.wvz[n];
                   TM.vvrSum += TM.mInertX[0]*wvBx*wvBx + TM.mInertY[0]*wvBy*wvBy + TM.mInertZ[0]*wvBz*wvBz;
               }
               TM.rotationalTemperature = (TM.vvrSum/TM.nMol)/3.0;
               TM.rotKinEnergyVal = 0.5 * TM.vvrSum/TM.nMol;
              
        }else if(TM.ThermostatType == 3){
               vFac = Math.sqrt(1+(TM.tau*(((TM.targetTemperature)/TM.instantTemperature)-1)));
               TM.vvSum = 0.0;
               for(int n = 0; n < TM.nMol; n++){
                   TM.rvx[n] *= vFac;     TM.rvy[n] *= vFac;    TM.rvz[n] *= vFac;
                   TM.vvSum += TM.rvx[n]*TM.rvx[n] + TM.rvy[n]*TM.rvy[n] + TM.rvz[n]*TM.rvz[n];
               }
               TM.translationalTemperature = (TM.vvSum/TM.nMol)/3.0;
               TM.trzKinEnergyVal = 0.5 * TM.vvSum/TM.nMol;
               TM.vvrSum = 0.0;
               for(int n = 0; n < TM.nMol; n++){
                   TM.wvx[n] *= vFac;     TM.wvy[n] *= vFac;    TM.wvz[n] *= vFac;
                   wvBx = TM.rMatT[(n*9)+0]*TM.wvx[n] + TM.rMatT[(n*9)+1]*TM.wvy[n] + TM.rMatT[(n*9)+2]*TM.wvz[n];
                   wvBy = TM.rMatT[(n*9)+3]*TM.wvx[n] + TM.rMatT[(n*9)+4]*TM.wvy[n] + TM.rMatT[(n*9)+5]*TM.wvz[n];
                   wvBz = TM.rMatT[(n*9)+6]*TM.wvx[n] + TM.rMatT[(n*9)+7]*TM.wvy[n] + TM.rMatT[(n*9)+8]*TM.wvz[n];
                   TM.vvrSum += TM.mInertX[0]*wvBx*wvBx + TM.mInertY[0]*wvBy*wvBy + TM.mInertZ[0]*wvBz*wvBz;
               }
               TM.rotationalTemperature = (TM.vvrSum/TM.nMol)/3.0;
               TM.rotKinEnergyVal = 0.5 * TM.vvrSum/TM.nMol;
        }else if(TM.ThermostatType == 0){
              // The empyt set thermostat. No thermostat is a subset of all thermostats. Nature takes care of itself.
              TM.vvSum = 0.0;
              for(int n = 0; n < TM.nMol; n++){
                  TM.vvSum += TM.rvx[n]*TM.rvx[n] + TM.rvy[n]*TM.rvy[n] + TM.rvz[n]*TM.rvz[n];
              }
              TM.translationalTemperature = (TM.vvSum/TM.nMol)/3.0;
              TM.trzKinEnergyVal = 0.5 * TM.vvSum/TM.nMol;

              TM.vvrSum = 0.0;
              for(int n = 0; n < TM.nMol; n++){
                  wvBx = TM.rMatT[(n*9)+0]*TM.wvx[n] + TM.rMatT[(n*9)+1]*TM.wvy[n] + TM.rMatT[(n*9)+2]*TM.wvz[n];
                  wvBy = TM.rMatT[(n*9)+3]*TM.wvx[n] + TM.rMatT[(n*9)+4]*TM.wvy[n] + TM.rMatT[(n*9)+5]*TM.wvz[n];
                  wvBz = TM.rMatT[(n*9)+6]*TM.wvx[n] + TM.rMatT[(n*9)+7]*TM.wvy[n] + TM.rMatT[(n*9)+8]*TM.wvz[n];
                  TM.vvrSum += TM.mInertX[0]*wvBx*wvBx + TM.mInertY[0]*wvBy*wvBy + TM.mInertZ[0]*wvBz*wvBz;
              }
              TM.rotationalTemperature = (TM.vvrSum/TM.nMol)/3.0;
              TM.rotKinEnergyVal = 0.5 * TM.vvrSum/TM.nMol;              
        }else{
              TM.vvSum = 0.0;
              for(int n = 0; n < TM.nMol; n++){
                  TM.rvx[n] *= TM.rvxf[n];     TM.rvy[n] *= TM.rvyf[n];    TM.rvz[n] *= TM.rvzf[n];
                  TM.vvSum += TM.rvx[n]*TM.rvx[n] + TM.rvy[n]*TM.rvy[n] + TM.rvz[n]*TM.rvz[n];
              }
              TM.translationalTemperature = (TM.vvSum/TM.nMol)/3.0;
              TM.trzKinEnergyVal = 0.5 * TM.vvSum/TM.nMol;

              TM.vvrSum = 0.0;
              for(int n = 0; n < TM.nMol; n++){
                  TM.wvx[n] *= TM.rrvxf[n];     TM.wvy[n] *= TM.rrvyf[n];    TM.wvz[n] *= TM.rrvzf[n];
                  wvBx = TM.rMatT[(n*9)+0]*TM.wvx[n] + TM.rMatT[(n*9)+1]*TM.wvy[n] + TM.rMatT[(n*9)+2]*TM.wvz[n];
                  wvBy = TM.rMatT[(n*9)+3]*TM.wvx[n] + TM.rMatT[(n*9)+4]*TM.wvy[n] + TM.rMatT[(n*9)+5]*TM.wvz[n];
                  wvBz = TM.rMatT[(n*9)+6]*TM.wvx[n] + TM.rMatT[(n*9)+7]*TM.wvy[n] + TM.rMatT[(n*9)+8]*TM.wvz[n];
                  TM.vvrSum += TM.mInertX[0]*wvBx*wvBx + TM.mInertY[0]*wvBy*wvBy + TM.mInertZ[0]*wvBz*wvBz;
              }
              TM.rotationalTemperature = (TM.vvrSum/TM.nMol)/3.0;
              TM.rotKinEnergyVal = 0.5 * TM.vvrSum/TM.nMol;
        }
        TM.kinEnergyVal = TM.rotKinEnergyVal + TM.trzKinEnergyVal;
        TM.instantTemperature = TM.translationalTemperature + TM.rotationalTemperature;
        TM.totEnergyVal = TM.kinEnergyVal + TM.potEnergyVal;
    }
}
