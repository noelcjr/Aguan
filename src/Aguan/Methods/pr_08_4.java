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
public class pr_08_4{ //extends MD{ In this case we do not want to extend MD for fidelity in testing.
                      //All the functions will be present here as in pr_08_4.c in a way that reflects
                      //that code as close as posible.
    private TheMatrix TM;   public restart RR;  private input IN;   public double timeNow;
    private MD md;
    public pr_08_4(String[] args){
        TM = new TheMatrix();
        IN = new input(args[1]);
        md = new MD();
        IN.GetNameList_pr_08_4(TM);
        RR = new restart(args[2]);
        //System.out.println("check 1");
        TM.restartIn = args[2];
        RR.ReadRestartHeader(TM);
        TM.initMatrix();
        //System.out.println("check 2");
        RR.ReadRestartFile(TM);
        //System.out.println("check 3");
        //IN.GetNameList_pr_08_4(TM);
        SetParams();
        SetupJob();
        TM.moreCycles = true;
        while(TM.moreCycles){
             SingleStep(TM,RR);
             if(TM.stepCount >= TM.stepLimit){
                TM.moreCycles = false;
             }
        }
    }
    public void SingleStep(TheMatrix TM, restart RR){
        TM.stepCount++;
        timeNow = TM.stepCount * TM.deltaT;
        System.out.printf("stepCount=%d  timeNow=%f \n",TM.stepCount,timeNow);
        LeapfrogStep(1, TM);
        md.GenSiteCoord(TM); 
        md.ComputeSiteForces(TM);
        ComputeTorqs(TM);
        ApplyThermostat(TM);
        LeapfrogStep(2, TM);
        md.ApplyBoundaryCond(TM);
        md.EvalProps(TM);
        if(TM.stepCount % TM.stepAdjustTemp == 0 || TM.stepCount < TM.stepEquil && TM.stepCount % 100 == 0){
           AdjustTemp(TM);}
        md.AccumProps(TM, 1);
        if(TM.stepCount % TM.stepAvg == 0){
            md.AccumProps(TM, 2);
            PrintSummary(TM);
            md.AccumProps(TM, 0);
        }
    }
    public void SetupJob(){
        DefineMol(TM);
    }
    public void SetParams(){
        TM.sitesMol  = 4;
        TM.regionX = (1/Math.pow(TM.density,(1/3)))*TM.initUcellX; //VSCopy(v2, s1, v1) = (v2).x = (s1) * (v1).x
        TM.regionY = (1/Math.pow(TM.density,(1/3)))*TM.initUcellY;
        TM.regionZ = (1/Math.pow(TM.density,(1/3)))*TM.initUcellZ;
        TM.nMol = TM.initUcellX * TM.initUcellY * TM.initUcellZ; //nMol = VProd (initUcell) = x * y * z
        TM.velMag = Math.sqrt(TM.NDIM*(1-(1/TM.nMol))*TM.targetTemperature);//velMag=sqrt(NDIM*(1-1/nMol)*temperature)
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
        for(int n = 0; n < TM.nMol; n++){
            TM.rax[n] = TM.ray[n] = TM.raz[n] = 0.0;
            torqSx = torqSy = torqSz = 0.0;
            for(int j = 0; j < TM.sitesMol; j++){
                  TM.rax[n] = TM.rax[n] + TM.fxs[n*TM.sitesMol + j];     // Right assumption?  Translational accelerations are summed,
                  TM.ray[n] = TM.ray[n] + TM.fys[n*TM.sitesMol + j];     // Rotational accelerations are then obtained from this
                  TM.raz[n] = TM.raz[n] + TM.fzs[n*TM.sitesMol + j];     // translational accelerations. Shouldn't accelerations
                  drx = TM.rxs[n*TM.sitesMol + j] - TM.rx[n];            // be distributed between tranlation and rotation? Instead,
                  dry = TM.rys[n*TM.sitesMol + j] - TM.ry[n];            // what seems to occur is that all the accelerations are added
                  drz = TM.rzs[n*TM.sitesMol + j] - TM.rz[n];            // to translation and then a fraction of those are used for rotational
                  tx = dry*TM.fzs[n*TM.sitesMol + j] - drz*TM.fys[n*TM.sitesMol + j];          // accelerations. There seems to be a double count.
                  ty = drz*TM.fxs[n*TM.sitesMol + j] - drx*TM.fzs[n*TM.sitesMol + j];          // some of the acceleration that goes to rotation
                  tz = drx*TM.fys[n*TM.sitesMol + j] - dry*TM.fxs[n*TM.sitesMol + j];          // is already going to translation. Also, the rotational
                  torqSx =  torqSx + tx;     torqSy =  torqSy + ty;     torqSz =  torqSz + tz; // inertia tensor is calculated based on masses only
            }                                                                                  // not charges. is this right?  could this be fixed?
            waBx = TM.rMatT[(n*9)+0]*torqSx + TM.rMatT[(n*9)+1]*torqSy + TM.rMatT[(n*9)+2]*torqSz;
            waBy = TM.rMatT[(n*9)+3]*torqSx + TM.rMatT[(n*9)+4]*torqSy + TM.rMatT[(n*9)+5]*torqSz;
            waBz = TM.rMatT[(n*9)+6]*torqSx + TM.rMatT[(n*9)+7]*torqSy + TM.rMatT[(n*9)+8]*torqSz;
            waBx /= TM.mInertX;
            waBy /= TM.mInertY;
            waBz /= TM.mInertZ;
            TM.wax[n] = TM.rMatT[(n*9)+0]*waBx + TM.rMatT[(n*9)+3]*waBy + TM.rMatT[(n*9)+6]*waBz;
            TM.way[n] = TM.rMatT[(n*9)+1]*waBx + TM.rMatT[(n*9)+4]*waBy + TM.rMatT[(n*9)+7]*waBz;
            TM.waz[n] = TM.rMatT[(n*9)+2]*waBx + TM.rMatT[(n*9)+5]*waBy + TM.rMatT[(n*9)+8]*waBz;
        }
    }
    public void DefineMol(TheMatrix TM){
       for(int j = 0; j < TM.sitesMol; j++){
          TM.rmx[j] = 0;   TM.rmy[j] = 0;  TM.rmz[0] = 0;
       }
       TM.rmz[0] =-0.0206;
       TM.rmz[1] = 0.0274;
       TM.rmy[2] = 0.240;    
       TM.rmz[2] = 0.165;
       TM.rmy[3] = -0.240;   
       TM.rmz[3] = 0.165;
       TM.mInertX = 0.00980;
       TM.mInertY = 0.00340;
       TM.mInertZ = 0.00640;
       TM.bCon = 183.5;
       TM.typeF[0] = 1;   
       TM.typeF[1] = 2;    
       TM.typeF[2] = 3;    
       TM.typeF[3] = 3;
    }
    public void ApplyThermostat (TheMatrix TM){
        double[] mc, mt;  //RMat mc, mt;
        mc = new double[9];      mt = new double[9];
        double vtx, vty, vtz, waBx, waBy, waBz, wvBx, wvBy, wvBz; //VecR vt, waB, wvB;
        double s1, s2, vFac;  //real s1, s2, vFac;
        int n;                //int n;

        s1 = s2 = 0;          //s1 = s2 = 0.;
        for(n = 0; n < TM.nMol; n++){         //DO_MOL {
        //VSAdd (vt, mol[n].rv, 0.5 * deltaT, mol[n].ra);
        vtx = TM.rvx[n] + (0.5*TM.deltaT*TM.rax[n]);
        vty = TM.rvy[n] + (0.5*TM.deltaT*TM.ray[n]);
        vtz = TM.rvz[n] + (0.5*TM.deltaT*TM.raz[n]);
        //s1 += VDot (vt, mol[n].ra);
        s1 += (vtx*TM.rax[n])+(vty*TM.ray[n])*(vtz*TM.raz[n]);
        //s2 += VLenSq (vt);
        s2 +=  (vtx*vtx)+(vty*vty)+(vtz*vtz);
        //VSAdd (vt, mol[n].wv, 0.5 * deltaT, mol[n].wa);
        vtx = TM.wvx[n] + (0.5*TM.deltaT*TM.wax[n]);
        vty = TM.wvy[n] + (0.5*TM.deltaT*TM.way[n]);
        vtz = TM.wvz[n] + (0.5*TM.deltaT*TM.waz[n]); 
        //MVMulT (wvB, mol[n].rMatT.u, vt);
        wvBx = (TM.rMatT[(n*9)+0]*vtx)+(TM.rMatT[(n*9)+1]*vty)+(TM.rMatT[(n*9)+2]*vtz);
        wvBy = (TM.rMatT[(n*9)+3]*vtx)+(TM.rMatT[(n*9)+4]*vty)+(TM.rMatT[(n*9)+5]*vtz);
        wvBz = (TM.rMatT[(n*9)+6]*vtx)+(TM.rMatT[(n*9)+7]*vty)+(TM.rMatT[(n*9)+8]*vtz);
        //MVMulT (waB, mol[n].rMatT.u, mol[n].wa);
        waBx = (TM.rMatT[(n*9)+0]*TM.wax[n])+(TM.rMatT[(n*9)+1]*TM.way[n])+(TM.rMatT[(n*9)+2]*TM.waz[n]);
        waBy = (TM.rMatT[(n*9)+3]*TM.wax[n])+(TM.rMatT[(n*9)+4]*TM.way[n])+(TM.rMatT[(n*9)+5]*TM.waz[n]);
        waBz = (TM.rMatT[(n*9)+6]*TM.wax[n])+(TM.rMatT[(n*9)+7]*TM.way[n])+(TM.rMatT[(n*9)+8]*TM.waz[n]); 
        //s1 += VWDot (mInert, wvB, waB);
        s1 += (TM.mInertX*wvBx*waBx)+(TM.mInertY*wvBy*waBy)+(TM.mInertZ*wvBz*waBz);
        //s2 += VWLenSq (mInert, wvB);
        s2 += (TM.mInertX*wvBx*wvBx)+(TM.mInertY*wvBy*wvBy)+(TM.mInertZ*wvBz*wvBz);
        } //}
        vFac = -s1 / s2;     //vFac = - s1 / s2;
        for(n = 0; n <TM.nMol; n++){  // DO_MOL {
        //VSAdd (vt, mol[n].rv, 0.5 * deltaT, mol[n].ra);
        vtx = TM.rvx[n] + (0.5*TM.deltaT*TM.rax[n]);
        vty = TM.rvy[n] + (0.5*TM.deltaT*TM.ray[n]);
        vtz = TM.rvz[n] + (0.5*TM.deltaT*TM.raz[n]);
        //VVSAdd (mol[n].ra, vFac, vt);
        TM.rax[n] = TM.rax[n] + (vFac*vtx);
        TM.ray[n] = TM.ray[n] + (vFac*vty);
        TM.raz[n] = TM.raz[n] + (vFac*vtz);
        //VSAdd (vt, mol[n].wv, 0.5 * deltaT, mol[n].wa);
        vtx = TM.wvx[n] + (0.5*TM.deltaT*TM.wax[n]);
        vty = TM.wvy[n] + (0.5*TM.deltaT*TM.way[n]);
        vtz = TM.wvz[n] + (0.5*TM.deltaT*TM.waz[n]);
        //VVSAdd (mol[n].wa, vFac, vt);
        TM.wax[n] = TM.wax[n] + (vFac*vtx);
        TM.way[n] = TM.way[n] + (vFac*vty);
        TM.waz[n] = TM.waz[n] + (vFac*vtz);
        }  // }
    }
    public void AdjustTemp(TheMatrix TM){
        double wvBx, wvBy, wvBz;
        double vFac;
        int n;
        
        TM.vvSum = 0.0;
        // DO_MOL vvSum += VLenSq (mol[n].rv);
        for(n = 0; n < TM.nMol; n++){
            TM.vvSum += TM.rvx[n]*TM.rvx[n] + TM.rvy[n]*TM.rvy[n] + TM.rvz[n]*TM.rvz[n];
        }
        // vFac = velMag / sqrt (vvSum / nMol);
        vFac = TM.velMag/Math.sqrt(TM.vvSum/TM.nMol);
        // DO_MOL VScale (mol[n].rv, vFac)
        for(n = 0; n < TM.nMol; n++){
            TM.rvx[n] *= vFac;  TM.rvy[n] *= vFac;   TM.rvz[n] *= vFac;
        }
        // vvSum = 0;
        TM.vvSum = 0.0;
        for(n = 0; n < TM.nMol; n++){
            // MVMulT (wvB, mol[n].rMatT.u, mol[n].wv);
            wvBx = TM.rMatT[(n*9)+0]*TM.wvx[n] + TM.rMatT[(n*9)+1]*TM.wvy[n] + TM.rMatT[(n*9)+2]*TM.wvz[n];
            wvBy = TM.rMatT[(n*9)+3]*TM.wvx[n] + TM.rMatT[(n*9)+4]*TM.wvy[n] + TM.rMatT[(n*9)+5]*TM.wvz[n];
            wvBz = TM.rMatT[(n*9)+6]*TM.wvx[n] + TM.rMatT[(n*9)+7]*TM.wvy[n] + TM.rMatT[(n*9)+8]*TM.wvz[n];
            // vvSum += VWLenSq (mInert, wvB);
            TM.vvSum += TM.mInertX*wvBx*wvBx + TM.mInertY*wvBy*wvBy + TM.mInertZ*wvBz*wvBz;
        }
        vFac = TM.velMag/Math.sqrt(TM.vvSum/TM.nMol);
        for(n = 0; n < TM.nMol; n++){
            TM.wvx[n] *= vFac;     TM.wvy[n] *= vFac;    TM.wvz[n] *= vFac;
        }
    }
    void PrintSummary(TheMatrix TM){
        System.out.printf("%5d %8.4f ",TM.stepCount, timeNow);
        System.out.printf("%7.4f ",(TM.vSumX+TM.vSumY+TM.vSumZ)/TM.nMol);
        System.out.printf("%7.4f %7.4f ",TM.potEnergyVal,TM.potEnergySum);
        System.out.printf("%7.4f %7.4f\n",TM.kinEnergyVal, TM.kinEnergySum);
    }
}
