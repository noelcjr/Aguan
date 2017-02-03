/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package Aguan.Methods;
import Aguan.TheMatrix.TheMatrix;
import Aguan.files.*;
import java.text.*;
/**
 *
 * @author Noel
 */
public class rotateTranslateOneWaterTest {
    private crd CRD;
    private psf PSF;
    private log LOG;
    public int regionX, regionY, regionZ;
    public rotateTranslateOneWaterTest(TheMatrix TM){
           CRD= new crd(TM.crd);
           PSF= new psf(TM.psf);
           LOG= new log(TM.logfile); 
           regionX = regionY = regionZ = 8;
           SetParams(TM);
           System.out.println("2 TM.nMol = "+TM.nMol);
           SetupJob(TM);
           System.out.println("3 TM.nMol = "+TM.nMol);
           // TIP3P only
           TM.rmz[0] = -0.0207;
           TM.rmz[1] = -0.0207;
           TM.rmy[2] = 0.2402;   TM.rmz[2] = 0.1653;
           TM.rmy[3] =-0.2402;   TM.rmz[3] = 0.1653;
    }
    public void SetParams(TheMatrix TM){
        if(TM.molType.equalsIgnoreCase("tip3p")){
           TM.sitesMol[0] = 4;
        }else if(TM.molType.equalsIgnoreCase("tip4p")){
           TM.sitesMol[0] = 4;
        }else if(TM.molType.equalsIgnoreCase("tip5p")){
           TM.sitesMol[0] = 5;
        }else if(TM.molType.equalsIgnoreCase("st2")){
           TM.sitesMol[0] = 5;
        }
        TM.NDIM = 3;
        TM.regionX = (1/Math.pow(TM.density,0.333333))*regionX;
        TM.regionY = (1/Math.pow(TM.density,0.333333))*regionY;
        TM.regionZ = (1/Math.pow(TM.density,0.333333))*regionZ;
 //     TM.nMol = regionX*regionY*regionZ;
        LOG.lout.println("Box = "+TM.regionX+" "+TM.regionY+" "+TM.regionZ);
    }
    public void SetupJob(TheMatrix TM){
        System.out.println("TEST 1: genereating molecule.");
        TM.initMatrix();
        DefineMol(TM,TM.molType);
        RandomGenerate(TM);
        GenSiteCoord(TM);
 //     TM.rvx[0] = 1.0;
        //PM.dcdM.write_header("CORD",(TM.minimization+TM.stepLimit),1, 1, (TM.nCharges[0]*TM.nMol));
        System.out.println("10 ="+TM.rx[0]+" "+TM.ry[0]+" "+TM.rz[0]);
        System.out.println("1w1 ="+TM.rxs[0]+" "+TM.rys[0]+" "+TM.rzs[0]);
        System.out.println("1w2 ="+TM.rxs[1]+" "+TM.rys[1]+" "+TM.rzs[1]);
        System.out.println("1w3 ="+TM.rxs[2]+" "+TM.rys[2]+" "+TM.rzs[2]);
        System.out.println("1w4 ="+TM.rxs[3]+" "+TM.rys[3]+" "+TM.rzs[3]);
/**     for(int i = 0; i < 2000; i++){
              TM.rx[0] = TM.rx[0] + TM.deltaT * TM.rvx[0];
              TM.ry[0] = TM.ry[0] + TM.deltaT * TM.rvy[0];
              TM.rz[0] = TM.rz[0] + TM.deltaT * TM.rvz[0];
              if(i < 10){System.out.println((i+1)+" ="+TM.rx[0]+" "+TM.ry[0]+" "+TM.rz[0]);}
              ApplyBoundaryCond(TM);
              GenSiteCoord(TM);
              PM.dcdM.write_dcdStep((TM.nCharges[0]*TM.nMol),TM.rxs,TM.rys,TM.rzs,TM.ro[0],TM.vdwPoint);
        }*/
        double[] mc, mt;
        double tx, ty, tz;
        mc = new double[9];          mt = new double[9];
        TM.wvx[0] = 0.0;  TM.wvy[0] = 0.0; TM.wvz[0] = 0.1;
        for(int n = 0; n < 200; n++){
            tx = TM.wvx[0];
            ty = TM.wvy[0];
            tz = TM.wvz[0];
            System.out.println("tx,ty,tz="+tx+" "+ty+" "+tz+"   Delta="+TM.deltaT);
            System.out.println(mc[0]+" "+mc[3]+" "+mc[6]);
            System.out.println(mc[1]+" "+mc[4]+" "+mc[7]);
            System.out.println(mc[2]+" "+mc[5]+" "+mc[8]);
            BuildStepRmat(mc, tx, ty, tz);
            System.out.println(" |"+mc[0]+" "+mc[3]+" "+mc[6]);
            System.out.println(" |"+mc[1]+" "+mc[4]+" "+mc[7]);
            System.out.println(" |"+mc[2]+" "+mc[5]+" "+mc[8]);
            MulMat(mt, mc, TM.rMatT, 3, 0);
            System.out.println("|"+mt[0]+" "+mt[3]+" "+mt[6]+"|");
            System.out.println("|"+mt[1]+" "+mt[4]+" "+mt[7]+"|");
            System.out.println("|"+mt[2]+" "+mt[5]+" "+mt[8]+"|");
            System.out.println("=");
            System.out.println("|"+mc[0]+" "+mc[3]+" "+mc[6]+"|");
            System.out.println("|"+mc[1]+" "+mc[4]+" "+mc[7]+"|");
            System.out.println("|"+mc[2]+" "+mc[5]+" "+mc[8]+"|");
            System.out.println("*");
            System.out.println("|"+TM.rMatT[0]+" "+TM.rMatT[3]+" "+TM.rMatT[6]+"|");
            System.out.println("|"+TM.rMatT[1]+" "+TM.rMatT[4]+" "+TM.rMatT[7]+"|");
            System.out.println("|"+TM.rMatT[2]+" "+TM.rMatT[5]+" "+TM.rMatT[8]+"|");
            TM.rMatT[0] = mt[0];    TM.rMatT[1] = mt[1];    TM.rMatT[2] = mt[2];
            TM.rMatT[3] = mt[3];    TM.rMatT[4] = mt[4];    TM.rMatT[5] = mt[5];
            TM.rMatT[6] = mt[6];    TM.rMatT[7] = mt[7];    TM.rMatT[8] = mt[8];
            ApplyBoundaryCond(TM);
            GenSiteCoord(TM);
            System.out.println("ww1 ="+TM.rxs[0]+" "+TM.rys[0]+" "+TM.rzs[0]);
            System.out.println("ww2 ="+TM.rxs[1]+" "+TM.rys[1]+" "+TM.rzs[1]);
            System.out.println("ww3 ="+TM.rxs[2]+" "+TM.rys[2]+" "+TM.rzs[2]);
            System.out.println("ww4 ="+TM.rxs[3]+" "+TM.rys[3]+" "+TM.rzs[3]);
            //PM.dcdM.write_dcdStep((TM.nCharges[0]*TM.nMol),TM.rxs,TM.rys,TM.rzs,TM.ro[0],TM.vdwPoint);
        }
        CRD.createCRDFile(TM);
        PSF.createPSFFile(TM);
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
            System.out.println("      c["+k+"] = "+c[k]+" s["+k+"]="+s[k]);
        }
        c0c2 = c[0] * c[2];
        c0s2 = c[0] * s[2];
        s0c2 = s[0] * c[2];
        s0s2 = s[0] * s[2];
        System.out.println("          c0c2="+c0c2+"  c0s2="+c0s2+"  s0c2="+s0c2+"  s0s2="+s0s2);
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
    public void ApplyBoundaryCond(TheMatrix TM){
      for(int n = 0; n < TM.nMol; n++){
           if( TM.rx[n] >= 0.5 * TM.regionX){TM.rx[n] = TM.rx[n] - TM.regionX;}
           else if(TM.rx[n] < -0.5 * TM.regionX){TM.rx[n] = TM.rx[n] + TM.regionX;}
           if(TM.ry[n] >= 0.5 * TM.regionY){TM.ry[n]= TM.ry[n]- TM.regionY;}
           else if(TM.ry[n] < -0.5 * TM.regionY){TM.ry[n] = TM.ry[n] + TM.regionY;}
           if(TM.rz[n] >= 0.5 * TM.regionZ){TM.rz[n] = TM.rz[n] - TM.regionZ;}
           else if(TM.rz[n] < -0.5 * TM.regionZ){TM.rz[n] = TM.rz[n] + TM.regionZ;}
       }
    }
    public void DefineMol(TheMatrix TM, String name){
        for(int j = 0; j < TM.sitesMol[0]; j++){
            TM.rmx[j] = 0.0; TM.rmy[j] = 0.0; TM.rmz[j] = 0.0;
        }
        if(name.equalsIgnoreCase("tip3p")){
           TM.bCon[0] = 120.3995;
           TM.ro[0] = 3.15061;
           TM.ep[0] = 0.1521;
           TM.rmz[0] = -0.0207;
           TM.rmz[1] = -0.0207;
           TM.rmy[2] = 0.2402;   TM.rmz[2] = 0.1653;
           TM.rmy[3] =-0.2402;   TM.rmz[3] = 0.1653;
           TM.typeF[0] = 1;   TM.typeF[1] = 2;    TM.typeF[2] = 3;    TM.typeF[3] = 3;
           TM.nAtoms[0] = 3;     TM.nCharges[0] = 3;
           // Only for tip4p it would have to be modified for other molecules or combination of molecules
           for(int j = 0; j < TM.nMol*TM.sitesMol[0]; j=j+4){
               TM.atomType[j] = 1;
               TM.atomType[j+1] = 2;
               TM.atomType[j+2] = 3;
               TM.atomType[j+3] = 3;
           }
           for(int j = 0; j < TM.nMol; j++){TM.moleculeIndex[j] = (TM.sitesMol[0]*3*j);}
        }else if(name.equalsIgnoreCase("tip4p")){
           TM.bCon[0] = 183.6615;
           TM.ro[0] = 3.15365;
           TM.ep[0] = 0.1549;
           TM.rmz[0] =-0.02065;
           TM.rmz[1] = 0.02695;
           TM.rmy[2] = 0.2400;    TM.rmz[2] = 0.1652;
           TM.rmy[3] = -0.2400;   TM.rmz[3] = 0.1652;
           TM.typeF[0] = 1;   TM.typeF[1] = 2;    TM.typeF[2] = 3;    TM.typeF[3] = 3;
           TM.nAtoms[0] = 3;      TM.nCharges[0] = 3;
           // Only for tip4p it would have to be modified for other molecules or combination of molecules
           for(int j = 0; j < TM.nMol*TM.sitesMol[0]; j=j+4){
               TM.atomType[j] = 1;
               TM.atomType[j+1] = 2;
               TM.atomType[j+2] = 3;
               TM.atomType[j+3] = 3;
           }
           for(int j = 0; j < TM.nMol; j++){TM.moleculeIndex[j] = (TM.sitesMol[0]*3*j);}
        }else if(name.equalsIgnoreCase("tip5p")){
           TM.bCon[0] = 39.4567;
           TM.ro[0] = 3.1200;
           TM.ep[0] = 0.1599;
           TM.rmz[0] = -0.0209;
           TM.rmy[1] = 0.2426;   TM.rmz[1] = 0.1669;
           TM.rmy[2] =-0.2426;   TM.rmz[2] = 0.1669;
           TM.rmx[3] = 0.1832;   TM.rmz[3] = -0.1504;
           TM.rmx[4] =-0.1832;   TM.rmz[4] = -0.1504;
           TM.typeF[0] = 1;   TM.typeF[1] = 3;    TM.typeF[2] = 3;
           TM.typeF[3] = 7;   TM.typeF[4] = 7;
           TM.nAtoms[0] = 3;     TM.nCharges[0] = 4;
           // Only for tip4p it would have to be modified for other molecules or combination of molecules
           for(int j = 0; j < TM.nMol*TM.sitesMol[0]; j=j+5){
               TM.atomType[j] = 2;
               TM.atomType[j+1] = 3;
               TM.atomType[j+2] = 3;
               TM.atomType[j+3] = 4;
               TM.atomType[j+4] = 4;
           }
           for(int j = 0; j < TM.nMol; j++){TM.moleculeIndex[j] = (TM.sitesMol[0]*3*j);}
        }else if(name.equalsIgnoreCase("st2")){
           TM.bCon[0] = 4.5224;
           TM.ro[0] = 3.1000;
           TM.ep[0] =0.0757;
           TM.rmx[0] = -0.0207;
           TM.rmy[1] = 0.263;   TM.rmz[1] = 0.1656;
           TM.rmy[2] =-0.263;   TM.rmz[2] = 0.1656;
           TM.rmx[3] = 0.2107;  TM.rmz[3] =-0.1697;
           TM.rmx[4] =-0.2107;  TM.rmz[4] =-0.1697;
           TM.typeF[0] = 1;   TM.typeF[1] = 3;    TM.typeF[2] = 3;    TM.typeF[3] = 7;   TM.typeF[4] = 7;
           TM.nAtoms[0] = 3;     TM.nCharges[0] = 4;
           // Only for tip4p it would have to be modified for other molecules or combination of molecules
           for(int j = 0; j < TM.nMol*TM.sitesMol[0]; j=j+5){
               TM.atomType[j] = 2;
               TM.atomType[j+1] = 3;
               TM.atomType[j+2] = 3;
               TM.atomType[j+3] = 4;
               TM.atomType[j+4] = 4;
           }
           for(int j = 0; j < TM.nMol; j++){TM.moleculeIndex[j] = (TM.sitesMol[0]*3*j);}
        }else if(name.equalsIgnoreCase("spc")){

        }
    }
    public void GenSiteCoord(TheMatrix TM){
        double tx, ty, tz;
        //System.out.println("5 TM.nMol = "+TM.nMol);
        for(int n = 0; n < TM.nMol; n++){
            for(int j = 0; j < TM.sitesMol[0]; j++){
             // System.out.println(n+" -- "+j);
                tx = TM.rMatT[(n*9)+0]*TM.rmx[j] + TM.rMatT[(n*9)+3]*TM.rmy[j] + TM.rMatT[(n*9)+6]*TM.rmz[j];
                ty = TM.rMatT[(n*9)+1]*TM.rmx[j] + TM.rMatT[(n*9)+4]*TM.rmy[j] + TM.rMatT[(n*9)+7]*TM.rmz[j];
                tz = TM.rMatT[(n*9)+2]*TM.rmx[j] + TM.rMatT[(n*9)+5]*TM.rmy[j] + TM.rMatT[(n*9)+8]*TM.rmz[j];

                TM.rxs[n*TM.sitesMol[0] + j] = TM.rx[n] + tx;
                TM.rys[n*TM.sitesMol[0] + j] = TM.ry[n] + ty;
                TM.rzs[n*TM.sitesMol[0] + j] = TM.rz[n] + tz;                
            }
        }
    }
    public void RandomGenerate(TheMatrix TM){
        double[] p, tq;
        double s = 0;
        int k, k1, k2;
        p = new double[10];  tq = new double[4];
        double a1, a2, a3, f;
        double eAngx, eAngy, eAngz;
        
        TM.rx[0] = 0.0;   TM.ry[0] = 0.0;    TM.rz[0] = 0.0;        
        TM.rvx[0] = 0.0;  TM.rvy[0] = 0.0;   TM.rvz[0] = 0.0;
        TM.rax[0] = 0.0;  TM.ray[0] = 0.0;   TM.raz[0] = 0.0;
        eAngx = 0.0;      eAngy = 0.0;       eAngz = 0.0;  
        a1 = 0.5 * eAngy;
        a2 = 0.5 * (eAngx - eAngz);
        a3 = 0.5 * (eAngx + eAngz);
        System.out.println("   a123 = "+a1+" "+a2+" "+a3);
        TM.q_u1[0] = Math.sin(a1) * Math.cos(a2);
        TM.q_u2[0] = Math.sin(a1) * Math.sin(a2);
        TM.q_u3[0] = Math.cos(a1) * Math.sin(a3);
        TM.q_u4[0] = Math.cos(a1) * Math.cos(a3);
        System.out.println("   q1234 = "+TM.q_u1[0]+" "+TM.q_u2[0]+" "+TM.q_u3[0]+" "+TM.q_u4[0]);
        tq[0] = TM.q_u1[0];  tq[1] = TM.q_u2[0];  tq[2] = TM.q_u3[0];    tq[3] = TM.q_u4[0];
        for(k = 0, k2 = 0; k2 < 4; k2++){
            for(k1 = k2; k1 < 4; k1++, k++){
                p[k] = 2.0*tq[k1]*tq[k2];
            }
        }
        TM.rMatT[0] = p[0] + p[9] - 1;  TM.rMatT[4] = p[4] + p[9] - 1;   TM.rMatT[8] = p[7] + p[9] - 1;
        s = 1.0;    //Transpose = 1
        TM.rMatT[1] = p[1] + s * p[8];   TM.rMatT[3] = p[1] - s * p[8];   TM.rMatT[2] = p[2] - s * p[6];
        TM.rMatT[6] = p[2] + s * p[6];   TM.rMatT[5] = p[5] + s * p[3];   TM.rMatT[7] = p[5] - s * p[3];
        System.out.println("    |"+TM.rMatT[0]+" "+TM.rMatT[3]+" "+TM.rMatT[6]+"|");
        System.out.println("    |"+TM.rMatT[1]+" "+TM.rMatT[4]+" "+TM.rMatT[7]+"|");
        System.out.println("    |"+TM.rMatT[2]+" "+TM.rMatT[5]+" "+TM.rMatT[8]+"|");
        
        TM.wvx[0] = 0.0;        TM.wvy[0] = 0.0;        TM.wvz[0] = 0.0;        
    }
}
