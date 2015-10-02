/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package Aguan.ntraj;

import Aguan.TheMatrix.TheMatrix;
import Aguan.files.*;
import java.text.*;
/**
 *
 * @author Noel
 */
public class generateMolecularOrientationTest {
    public int regionX, regionY, regionZ;
    private int counter = 0;
    private int c0, c1, c2, c3, c4, c5, c6, c7;
    private int tests, runs;
    private crd CRD;
    private psf PSF;
    private log LOG;
    public generateMolecularOrientationTest(TheMatrix TM){
           CRD = new crd(TM.crd); 
           PSF = new psf(TM.psf);   
           LOG = new log(TM.logfile);       
           c0 = c1 = c2 = c3 = c4 = c5 = c6 = c7 = 0;
           regionX = regionY = regionZ = 8;
           SetParams(TM);
           SetupJob(TM);
    }
    public void SetParams(TheMatrix TM){
        if(TM.molType.equalsIgnoreCase("tip3p")){
           TM.sitesMol = 4;
        }else if(TM.molType.equalsIgnoreCase("tip4p")){
           TM.sitesMol = 4;
        }else if(TM.molType.equalsIgnoreCase("tip5p")){
           TM.sitesMol = 5;
        }else if(TM.molType.equalsIgnoreCase("st2")){
           TM.sitesMol = 5;
        }
        TM.NDIM = 3;
        TM.regionX = (1/Math.pow(TM.density,0.333333))*regionX;
        TM.regionY = (1/Math.pow(TM.density,0.333333))*regionY;
        TM.regionZ = (1/Math.pow(TM.density,0.333333))*regionZ;
        TM.nMol = 8;
        LOG.lout.println("Box = "+TM.regionX+" "+TM.regionY+" "+TM.regionZ);
    }
    public void SetupJob(TheMatrix TM){
        System.out.println("TEST 1: genereating molecules in 45 degree ZXZ rotation differences.");
        TM.initMatrix();
        DefineMol(TM,TM.molType);
   //     GenerateAllXYZ45(TM,PM);
        GenerateAllXYZ45_2(TM);
        System.out.println("Moltype:"+TM.molType);
        System.out.println("nMol:"+TM.nMol);
        System.out.println("nMol:"+TM.sitesMol);
 //     testRotMatrix(TM,PM);
        CRD.createCRDFile(TM);
        PSF.createPSFFile(TM);
        tests = 1000;  runs = 100000;
  //    System.out.println("TEST 2: genereating and histograming Random vectors inton 10 bins.");
  //    generateRandomNumbers();
  //    System.out.println("TEST 3: genereating and histograming Random vectors inton 10 bins.");
  //    testQuadrantFromRandomVectors();
  //    System.out.println("TEST 4: genereating Random vectors constructing water molecules and.");
  //    System.out.println("testing the quadrant they get allocated, and their distribution. TIP3P only.");
  //    testWaterGeneratinQuadrantAllocation(TM,PM);
    }
    public void testRotMatrix(TheMatrix TM){
           double[] p, tq;
           double s = 0;
           int k, k1, k2;
           p = new double[10];  tq = new double[4];
           double rx = 0; double ry = 0; double rz = 0;
           double tx, ty, tz;
           String pattern = "##0.00000";
           DecimalFormat Format = new DecimalFormat(pattern);
           double eAngx = Math.atan2(rx,ry);     
           double eAngy = Math.acos(rz);             
           double temp = Math.random();
           double eAngz = 2*Math.PI*temp;  
           double a1 = 0.5*(eAngy);
           double a2 = 0.5*(eAngx-eAngz);
           double a3 = 0.5*(eAngx+eAngz);
           double[] mc = new double[9];          double[] mt = new double[9];
           for(int n = 0; n < 50; n++){
                TM.rx[n] = TM.ry[n] = TM.rz[n] = 0.0;
                TM.rvx[n] = TM.rvy[n] = TM.rvz[n] = 0.0;
                TM.wvx[n] = 0.0;  TM.wvy[n] = 0.0; TM.wvz[n] = 0.0;
                TM.rax[n] = TM.ray[n] = TM.raz[n] = 0.0;
                TM.q_u1[n] = Math.sin(0)*Math.cos(0);
                TM.q_u2[n] = Math.sin(0)*Math.sin(0);
                TM.q_u3[n] = Math.cos(0)*Math.sin(0);
                TM.q_u4[n] = Math.cos(0)*Math.cos(0);  
                tq[0] = TM.q_u1[n];  tq[1] = TM.q_u2[n];  tq[2] = TM.q_u3[n];    tq[3] = TM.q_u4[n];
                for(k = 0, k2 = 0; k2 < 4; k2++){
                    for(k1 = k2; k1 < 4; k1++, k++){
                        p[k] = 2.0*tq[k1]*tq[k2];
                        LOG.lout.println("   p["+k+"]"+k1+","+k2+" = "+Format.format(p[k]));
                    }
                }
                TM.rMatT[n*9+0] = p[0] + p[9] - 1;  TM.rMatT[n*9+4] = p[4] + p[9] - 1;   TM.rMatT[n*9+8] = p[7] + p[9] - 1;
                s = 1.0;    //Transpose = 1
                TM.rMatT[n*9+1] = p[1] + s * p[8];   TM.rMatT[n*9+3] = p[1] - s * p[8];   TM.rMatT[n*9+2] = p[2] - s * p[6];
                TM.rMatT[n*9+6] = p[2] + s * p[6];   TM.rMatT[n*9+5] = p[5] + s * p[3];   TM.rMatT[n*9+7] = p[5] - s * p[3];
                LOG.lout.println("|"+TM.rMatT[0]+" "+TM.rMatT[1]+" "+TM.rMatT[2]+"|");
                LOG.lout.println("|"+TM.rMatT[3]+" "+TM.rMatT[4]+" "+TM.rMatT[5]+"|");
                LOG.lout.println("|"+TM.rMatT[6]+" "+TM.rMatT[7]+" "+TM.rMatT[8]+"|");
                for(int j = 0; j < TM.sitesMol; j++){
                    tx = TM.rMatT[(n*9)+0]*TM.rmx[j] + TM.rMatT[(n*9)+3]*TM.rmy[j] + TM.rMatT[(n*9)+6]*TM.rmz[j];
                    ty = TM.rMatT[(n*9)+1]*TM.rmx[j] + TM.rMatT[(n*9)+4]*TM.rmy[j] + TM.rMatT[(n*9)+7]*TM.rmz[j];
                    tz = TM.rMatT[(n*9)+2]*TM.rmx[j] + TM.rMatT[(n*9)+5]*TM.rmy[j] + TM.rMatT[(n*9)+8]*TM.rmz[j];

                    TM.rxs[n*TM.sitesMol + j] = TM.rx[n] + tx;
                    TM.rys[n*TM.sitesMol + j] = TM.ry[n] + ty;
                    TM.rzs[n*TM.sitesMol + j] = TM.rz[n] + tz;
                  //  LOG.lout.println(j+" = "+(TM.rxs[n*TM.sitesMol + j]*TM.ro)+" "+(TM.rys[n*TM.sitesMol + j]*TM.ro)+" "+(TM.rzs[n*TM.sitesMol + j]*TM.ro));
                }           
                tx = TM.wvx[n] + (n*0.1);
                ty = TM.wvy[n];
                tz = TM.wvz[n];
                BuildStepRmat(mc, tx, ty, tz);
                MulMat(mt, mc, TM.rMatT, 3, 0);
                TM.rMatT[n*9] = mt[0];    TM.rMatT[n*9+1] = mt[1];    TM.rMatT[n*9+2] = mt[2];
                TM.rMatT[n*9+3] = mt[3];    TM.rMatT[n*9+4] = mt[4];    TM.rMatT[n*9+5] = mt[5];
                TM.rMatT[n*9+6] = mt[6];    TM.rMatT[n*9+7] = mt[7];    TM.rMatT[n*9+8] = mt[8];
                for(int j = 0; j < TM.sitesMol; j++){   // tip3p only
                    tx = TM.rMatT[(n*9)+0]*TM.rmx[j] + TM.rMatT[(n*9)+3]*TM.rmy[j] + TM.rMatT[(n*9)+6]*TM.rmz[j];
                    ty = TM.rMatT[(n*9)+1]*TM.rmx[j] + TM.rMatT[(n*9)+4]*TM.rmy[j] + TM.rMatT[(n*9)+7]*TM.rmz[j];
                    tz = TM.rMatT[(n*9)+2]*TM.rmx[j] + TM.rMatT[(n*9)+5]*TM.rmy[j] + TM.rMatT[(n*9)+8]*TM.rmz[j];
  
                    TM.rxs[n*TM.sitesMol + j] = TM.rx[n] + tx;
                    TM.rys[n*TM.sitesMol + j] = TM.ry[n] + ty;
                    TM.rzs[n*TM.sitesMol + j] = TM.rz[n] + tz;
                }
         }  
    }
    public void GenerateAllXYZ45_2(TheMatrix TM){
          double[] p, tq;
          double s = 0;
          int k, k1, k2;
          p = new double[10];  tq = new double[4];
          double a1, a2, a3, f;
          double eAngx, eAngy, eAngz;
          String pattern = "##0.00000";
          DecimalFormat Format = new DecimalFormat(pattern);
          int n;
          n = 0;
          int phi, theta, psi;
          phi = theta = psi = 0;  
          double tx, ty, tz;
          for(theta = 0; theta < 360; theta=theta+45){
              TM.rx[n] = 0; 
              TM.ry[n] = 0; 
              TM.rz[n] = 0; 
              TM.q_u1[n] = Math.sin(0.5*theta*Math.PI/180) * Math.cos(0.5*(phi-psi)*Math.PI/180);
              TM.q_u2[n] = Math.sin(0.5*theta*Math.PI/180) * Math.sin(0.5*(phi-psi)*Math.PI/180);
              TM.q_u3[n] = Math.cos(0.5*theta*Math.PI/180) * Math.sin(0.5*(phi+psi)*Math.PI/180);
              TM.q_u4[n] = Math.cos(0.5*theta*Math.PI/180) * Math.cos(0.5*(phi+psi)*Math.PI/180);
              tq[0] = TM.q_u1[n];  tq[1] = TM.q_u2[n];  tq[2] = TM.q_u3[n];    tq[3] = TM.q_u4[n];
              for(k = 0, k2 = 0; k2 < 4; k2++){
                  for(k1 = k2; k1 < 4; k1++, k++){
                      p[k] = 2.0*tq[k1]*tq[k2];
                  }
              }
              TM.rMatT[n*9+0] = p[0] + p[9] - 1;  TM.rMatT[n*9+4] = p[4] + p[9] - 1;   TM.rMatT[n*9+8] = p[7] + p[9] - 1;
              s = 1.0;    //Transpose = 1
              TM.rMatT[n*9+1] = p[1] + s * p[8];   TM.rMatT[n*9+3] = p[1] - s * p[8];   TM.rMatT[n*9+2] = p[2] - s * p[6];
              TM.rMatT[n*9+6] = p[2] + s * p[6];   TM.rMatT[n*9+5] = p[5] + s * p[3];   TM.rMatT[n*9+7] = p[5] - s * p[3];
              for(int j = 0; j < TM.sitesMol; j++){
                  tx = TM.rMatT[(n*9)+0]*TM.rmx[j] + TM.rMatT[(n*9)+3]*TM.rmy[j] + TM.rMatT[(n*9)+6]*TM.rmz[j];
                  ty = TM.rMatT[(n*9)+1]*TM.rmx[j] + TM.rMatT[(n*9)+4]*TM.rmy[j] + TM.rMatT[(n*9)+7]*TM.rmz[j];
                  tz = TM.rMatT[(n*9)+2]*TM.rmx[j] + TM.rMatT[(n*9)+5]*TM.rmy[j] + TM.rMatT[(n*9)+8]*TM.rmz[j];

                  TM.rxs[n*TM.sitesMol + j] = TM.rx[n] + tx;
                  TM.rys[n*TM.sitesMol + j] = TM.ry[n] + ty;
                  TM.rzs[n*TM.sitesMol + j] = TM.rz[n] + tz;
              }                                                                   
              n++;
          }
          System.out.println("n="+n);
    }
    public void GenerateAllXYZ45(TheMatrix TM){
          double[] p, tq;
          double s = 0;
          int k, k1, k2;
          p = new double[10];  tq = new double[4];
          double a1, a2, a3, f;
          double eAngx, eAngy, eAngz;
          String pattern = "##0.00000";
          DecimalFormat Format = new DecimalFormat(pattern);
          int a, b, c, n;
          a = b = c = 1;
          n = 0;
          double tx, ty, tz;
          double distHHx, distHHy, distHHz;
          for(int x1 = 0; x1 < 360; x1=x1+45){
              for(int z1 = 0; z1 < 360; z1=z1+45){
                  for(int x2 = 0; x2 < 360; x2=x2+45){
                      TM.rx[n] = 0; //a*(TM.regionX/regionX);  // or zero
                      TM.ry[n] = 0; //b*(TM.regionY/regionY);  // or zero
                      TM.rz[n] = 0; //c*(TM.regionZ/regionZ);  // or zero
                      if(n < 100){
                          LOG.lout.print((n+1) +")  "+x1+" "+z1+" "+x2+" - "+a+" "+b+" "+c);
                       //   LOG.lout.println("    x,y,z");
                       //   LOG.lout.println("    "+Format.format(TM.rx[n])+" "+Format.format(TM.ry[n])+" "+Format.format(TM.rz[n]));
                          TM.q_u1[n] = Math.sin(0.5*x2*Math.PI/180) * Math.cos(0.5*(x1-z1)*Math.PI/180);
                          TM.q_u2[n] = Math.sin(0.5*x2*Math.PI/180) * Math.sin(0.5*(x1-z1)*Math.PI/180);
                          TM.q_u3[n] = Math.cos(0.5*x2*Math.PI/180) * Math.sin(0.5*(x1+z1)*Math.PI/180);
                          TM.q_u4[n] = Math.cos(0.5*x2*Math.PI/180) * Math.cos(0.5*(x1+z1)*Math.PI/180);
                       //   LOG.lout.println("    q1,2,3,4");
                       //   LOG.lout.println("    "+Format.format(TM.q_u1[n])+" "+Format.format(TM.q_u2[n])+" "+Format.format(TM.q_u3[n])+" "+Format.format(TM.q_u4[n]));
                          tq[0] = TM.q_u1[n];  tq[1] = TM.q_u2[n];  tq[2] = TM.q_u3[n];    tq[3] = TM.q_u4[n];
                          for(k = 0, k2 = 0; k2 < 4; k2++){
                              for(k1 = k2; k1 < 4; k1++, k++){
                                  p[k] = 2.0*tq[k1]*tq[k2];
                       //           LOG.lout.println("   p["+k+"]"+k1+","+k2+" = "+Format.format(p[k]));
                              }
                          }
                          TM.rMatT[n*9+0] = p[0] + p[9] - 1;  TM.rMatT[n*9+4] = p[4] + p[9] - 1;   TM.rMatT[n*9+8] = p[7] + p[9] - 1;
                          s = 1.0;    //Transpose = 1
                          TM.rMatT[n*9+1] = p[1] + s * p[8];   TM.rMatT[n*9+3] = p[1] - s * p[8];   TM.rMatT[n*9+2] = p[2] - s * p[6];
                          TM.rMatT[n*9+6] = p[2] + s * p[6];   TM.rMatT[n*9+5] = p[5] + s * p[3];   TM.rMatT[n*9+7] = p[5] - s * p[3];
                          for(int j = 0; j < TM.sitesMol; j++){
                                tx = TM.rMatT[(n*9)+0]*TM.rmx[j] + TM.rMatT[(n*9)+3]*TM.rmy[j] + TM.rMatT[(n*9)+6]*TM.rmz[j];
                                ty = TM.rMatT[(n*9)+1]*TM.rmx[j] + TM.rMatT[(n*9)+4]*TM.rmy[j] + TM.rMatT[(n*9)+7]*TM.rmz[j];
                                tz = TM.rMatT[(n*9)+2]*TM.rmx[j] + TM.rMatT[(n*9)+5]*TM.rmy[j] + TM.rMatT[(n*9)+8]*TM.rmz[j];

                                TM.rxs[n*TM.sitesMol + j] = TM.rx[n] + tx;
                                TM.rys[n*TM.sitesMol + j] = TM.ry[n] + ty;
                                TM.rzs[n*TM.sitesMol + j] = TM.rz[n] + tz;
                          }
                          distHHx = (TM.rxs[n*TM.sitesMol+2]+TM.rxs[n*TM.sitesMol+3])/2;
                          distHHy = (TM.rys[n*TM.sitesMol+2]+TM.rys[n*TM.sitesMol+3])/2;
                          distHHz = (TM.rzs[n*TM.sitesMol+2]+TM.rzs[n*TM.sitesMol+3])/2;
                     //     LOG.lout.println("    (1,0,0)*|"+Format.format(TM.rMatT[n*9+0])+" "+Format.format(TM.rMatT[n*9+1])+" "+Format.format(TM.rMatT[n*9+2])+"|");
                     //     LOG.lout.println("            |"+Format.format(TM.rMatT[n*9+3])+" "+Format.format(TM.rMatT[n*9+4])+" "+Format.format(TM.rMatT[n*9+5])+"|");
                     //     LOG.lout.println("            |"+Format.format(TM.rMatT[n*9+6])+" "+Format.format(TM.rMatT[n*9+7])+" "+Format.format(TM.rMatT[n*9+8])+"|");
                     //     LOG.lout.println("    =("+Format.format((1*TM.rMatT[n*9+0]))+","+Format.format((1*TM.rMatT[n*9+1]))+","+Format.format((1*TM.rMatT[n*9+2]))+") => ");
                     //     LOG.lout.print("Quadrant = "+getQuadrant(1*TM.rMatT[n*9+0],1*TM.rMatT[n*9+1],1*TM.rMatT[n*9+2]));
                          LOG.lout.print("  Water Vector "+n+" = ("+Format.format(distHHx*TM.ro)+","+Format.format(distHHy*TM.ro)+","+Format.format(distHHz*TM.ro)+")");
                          LOG.lout.println("Q0 = "+getQuadrant(distHHx*TM.ro,distHHy*TM.ro,distHHz*TM.ro)+"   Q1 = "+getQuadrant1(distHHx*TM.ro,distHHy*TM.ro,distHHz*TM.ro));
               //           LOG.lout.println("   Quadrant2 = "+getQuadrant2(distHHx*TM.ro,distHHy*TM.ro,distHHz*TM.ro));
                          LOG.lout.println("-------------------------------------------------------------------------------------------------");
                      }else{
                          LOG.lout.print(n +") "+x1+" "+z1+" "+x2+" - "+a+" "+b+" "+c);
                          //  LOG.lout.println("    x,y,z");
                          //  LOG.lout.println("    "+Format.format(TM.rx[n])+" "+Format.format(TM.ry[n])+" "+Format.format(TM.rz[n]));
                          TM.q_u1[n] = Math.sin(0.5*x2*Math.PI/180) * Math.cos(0.5*(x1-z1)*Math.PI/180);
                          TM.q_u2[n] = Math.sin(0.5*x2*Math.PI/180) * Math.sin(0.5*(x1-z1)*Math.PI/180);
                          TM.q_u3[n] = Math.cos(0.5*x2*Math.PI/180) * Math.sin(0.5*(x1+z1)*Math.PI/180);
                          TM.q_u4[n] = Math.cos(0.5*x2*Math.PI/180) * Math.cos(0.5*(x1+z1)*Math.PI/180);
                          //  LOG.lout.println("    q1,2,3,4");
                          //  LOG.lout.println("    "+Format.format(TM.q_u1[n])+" "+Format.format(TM.q_u2[n])+" "+Format.format(TM.q_u3[n])+" "+Format.format(TM.q_u4[n]));

                          tq[0] = TM.q_u1[n];  tq[1] = TM.q_u2[n];  tq[2] = TM.q_u3[n];    tq[3] = TM.q_u4[n];
                          for(k = 0, k2 = 0; k2 < 4; k2++){
                              for(k1 = k2; k1 < 4; k1++, k++){
                                  p[k] = 2.0*tq[k1]*tq[k2];
                                 //   LOG.lout.println("   p["+k+"]"+k1+","+k2+" = "+p[k]);
                              }
                          }
                          TM.rMatT[n*9+0] = p[0] + p[9] - 1;  TM.rMatT[n*9+4] = p[4] + p[9] - 1;   TM.rMatT[n*9+8] = p[7] + p[9] - 1;
                          s = 1.0;    //Transpose = 1
                          TM.rMatT[n*9+1] = p[1] + s * p[8];   TM.rMatT[n*9+3] = p[1] - s * p[8];   TM.rMatT[n*9+2] = p[2] - s * p[6];
                          TM.rMatT[n*9+6] = p[2] + s * p[6];   TM.rMatT[n*9+5] = p[5] + s * p[3];   TM.rMatT[n*9+7] = p[5] - s * p[3];

                          for(int j = 0; j < TM.sitesMol; j++){
                                tx = TM.rMatT[(n*9)+0]*TM.rmx[j] + TM.rMatT[(n*9)+3]*TM.rmy[j] + TM.rMatT[(n*9)+6]*TM.rmz[j];
                                ty = TM.rMatT[(n*9)+1]*TM.rmx[j] + TM.rMatT[(n*9)+4]*TM.rmy[j] + TM.rMatT[(n*9)+7]*TM.rmz[j];
                                tz = TM.rMatT[(n*9)+2]*TM.rmx[j] + TM.rMatT[(n*9)+5]*TM.rmy[j] + TM.rMatT[(n*9)+8]*TM.rmz[j];

                                TM.rxs[n*TM.sitesMol + j] = TM.rx[n] + tx;
                                TM.rys[n*TM.sitesMol + j] = TM.ry[n] + ty;
                                TM.rzs[n*TM.sitesMol + j] = TM.rz[n] + tz;
                          }
                          distHHx = (TM.rxs[n*TM.sitesMol+2]+TM.rxs[n*TM.sitesMol+3])/2;
                          distHHy = (TM.rys[n*TM.sitesMol+2]+TM.rys[n*TM.sitesMol+3])/2;
                          distHHz = (TM.rzs[n*TM.sitesMol+2]+TM.rzs[n*TM.sitesMol+3])/2;
                          //     LOG.lout.println("    (1,0,0)*|"+Format.format(TM.rMatT[n*9+0])+" "+Format.format(TM.rMatT[n*9+1])+" "+Format.format(TM.rMatT[n*9+2])+"|");
                          //     LOG.lout.println("            |"+Format.format(TM.rMatT[n*9+3])+" "+Format.format(TM.rMatT[n*9+4])+" "+Format.format(TM.rMatT[n*9+5])+"|");
                          //     LOG.lout.println("            |"+Format.format(TM.rMatT[n*9+6])+" "+Format.format(TM.rMatT[n*9+7])+" "+Format.format(TM.rMatT[n*9+8])+"|");
                          //     LOG.lout.println("    =("+Format.format((1*TM.rMatT[n*9+0]))+","+Format.format((1*TM.rMatT[n*9+1]))+","+Format.format((1*TM.rMatT[n*9+2]))+") => ");
                          //     LOG.lout.print("Quadrant = "+getQuadrant(1*TM.rMatT[n*9+0],1*TM.rMatT[n*9+1],1*TM.rMatT[n*9+2]));
                          LOG.lout.print("  Water Vector "+n+" = ("+Format.format(distHHx*TM.ro)+","+Format.format(distHHy*TM.ro)+","+Format.format(distHHz*TM.ro)+")");
                          LOG.lout.println("Q0 = "+getQuadrant(distHHx*TM.ro,distHHy*TM.ro,distHHz*TM.ro)+"   Q1 = "+getQuadrant1(distHHx*TM.ro,distHHy*TM.ro,distHHz*TM.ro));
                     //   LOG.lout.println("   Quadrant2 = "+getQuadrant2(distHHx*TM.ro,distHHy*TM.ro,distHHz*TM.ro));
                          LOG.lout.println("-------------------------------------------------------------------------------------------------");
                      }
                      c++; n++;
                  }
                  b++; c = 1;
              }
              a++;
              b = 1;
          }
    }
    /**TEST: Random Number Generator. I wasn't sure how good it was or if
       I had to implement a better random number generator. The java standard random
       number generator seem to deviate 2, 3 or max 4 hundred per 1000000 vectors
       generated; that is, if a million vectors are generates and binned in ten
       groups evenly, the groups differed by 2 or 3 hundred vector per 100000.**/
    public void generateRandomNumbers(){
        double rx, ry, rz;        
        int[] ten = new int[10];
        float[] ave = new float[10];
        float[] aveDiff = new float[tests];
        ave[0] = ave[1] = ave[2] = ave[3] = ave[4] = ave[5] = ave[6] = ave[7] = ave[8] = ave[9] = 0;
        int total1, count;
        float total2;
        total2 = count = 0;
        float aveDiff_TOTAL = 0;
        int multi = 16384;
        for(int j = 0; j < tests; j++){
           total1 = 0;
           ten[0] = ten[1] = ten[2] = ten[3] = ten[4] = ten[5] = ten[6] = ten[7] = ten[8] = ten[9] = 0;
           for(int i = 0; i < multi; i++){
              rx = Math.random();
              if(rx > 0 && rx <= 0.1){ten[0]++;
              }else if(rx > 0.1 && rx <= 0.2){ten[1]++;
              }else if(rx > 0.2 && rx <= 0.3){ten[2]++;
              }else if(rx > 0.3 && rx <= 0.4){ten[3]++;
              }else if(rx > 0.4 && rx <= 0.5){ten[4]++;
              }else if(rx > 0.5 && rx <= 0.6){ten[5]++;
              }else if(rx > 0.6 && rx <= 0.7){ten[6]++;
              }else if(rx > 0.7 && rx <= 0.8){ten[7]++;
              }else if(rx > 0.8 && rx <= 0.9){ten[8]++;
              }else if(rx > 0.9 && rx <= 1.0){ten[9]++;}
           }
           ave[0] = ave[0] + ten[0];   ave[1] = ave[1] + ten[1];   ave[2] = ave[2] + ten[2];   ave[3] = ave[3] + ten[3];   ave[4] = ave[4] + ten[4];
           ave[5] = ave[5] + ten[5];   ave[6] = ave[6] + ten[6];   ave[7] = ave[7] + ten[7];   ave[8] = ave[8] + ten[8];   ave[9] = ave[9] + ten[9];
           total1 = ten[0] + ten[1] + ten[2] + ten[3] + ten[4] + ten[5] + ten[6] + ten[7] + ten[8] + ten[9];
           int countAVE = 1;
           float aveDiff_TEN = 0;
           for(int h=0;h < 10; h++){
               aveDiff_TEN = aveDiff_TEN + Math.abs((multi/10)-ten[h]);
               //  System.out.println(countAVE+" "+h+" "+hh+" "+ten[h]+" "+ten[hh]);
               countAVE++;               
           }
           aveDiff_TOTAL = aveDiff_TOTAL + (aveDiff_TEN/countAVE);
           System.out.print("j "+j+" = "+ten[0]+" "+ten[1]+" "+ten[2]+" "+ten[3]+" "+ten[4]+" "+ten[5]+" "+ten[6]+" "+ten[7]+" "+ten[8]+" "+ten[9]+" = "+total1);
           System.out.println("     Average bin difference = "+(aveDiff_TEN/countAVE)+" = "+((aveDiff_TEN/countAVE)/(multi)*100)+"%");
           count++;
        }
        ave[0] /= count; ave[1] /= count;  ave[2] /= count;  ave[3] /= count;  ave[4] /= count;
        ave[5] /= count; ave[6] /= count;  ave[7] /= count;  ave[8] /= count;  ave[9] /= count;
        total2 = ave[0] + ave[1] + ave[2] + ave[3] + ave[4] + ave[5] + ave[6] + ave[7] + ave[8] + ave[9];
        System.out.print("ave = "+ave[0]+" "+ave[1]+" "+ave[2]+" "+ave[3]+" "+ave[4]+" "+ave[5]+" "+ave[6]+" "+ave[7]+" "+ave[8]+" "+ave[9]+" = "+total2);
        System.out.println("     Total Average bin difference = "+(aveDiff_TOTAL/10)+"="+((aveDiff_TOTAL/10)/(multi)*100)+"%");
        System.out.println();
    }
    /**
     *
     * @param Rval
     * @param Rpl
     * @return
     * TEST: I generate random vectors and then I used getQuadrant2 to
     * figure out which quadrant the vector is in.
     */
    public void testQuadrantFromRandomVectors(){
        double rx, ry, rz;
        int total1;
        float total2;
        int q, q0, q1, q2, q3, q4, q5, q10, q21, q30, q32, q40, q41, q42, q43, q50, q51, q52, q53, q410, q421, q430, q432, q510, q521, q530, q532;
        q = q0 = q1 = q2 = q3 = q4 = q5 = 0;
        float[] aveQ = new float[26];
        aveQ[0] = aveQ[1] = aveQ[2] = aveQ[3] = aveQ[4] = aveQ[5] = 0;
        aveQ[6] = aveQ[7] = aveQ[8] = aveQ[9] = aveQ[10] = aveQ[11] = 0;
        aveQ[12] = aveQ[13] = aveQ[14] = aveQ[15] = aveQ[16] = aveQ[17] = 0;
        aveQ[18] = aveQ[19] = aveQ[20] = aveQ[21] = aveQ[22] = aveQ[23] = 0;
        aveQ[24] = aveQ[25] = 0;

        int count = 0;
        int multi = 16384;
        float aveDiff_TOTAL = 0;
        for(int k = 0; k < tests; k++){
           q0 = q1 = q2 = q3 = q4 = q5 = 0;
           q10 = q21 = q30 = q32 = q40 = 0;
           q41 = q42 = q43 = q50 = q51 = 0;
           q52 = q53 = q410 = q421 = q430 = 0;
           q432 = q510 = q521 = q530 = q532 = 0;
           for(int ll = 0; ll < multi; ll++){
              rx = Math.random();
              ry = Math.random();
              rz = Math.random();
              if(Math.random() < 0.5){rx=-1*rx;}
              if(Math.random() < 0.5){ry=-1*ry;}
              if(Math.random() < 0.5){rz=-1*rz;}              
              q = getQuadrant(rx,ry,rz);
                 if(q == 0){q0++;
              }else if(q == 1){q1++;
              }else if(q == 2){q2++;
              }else if(q == 3){q3++;
              }else if(q == 4){q4++;
              }else if(q == 5){q5++;
              }else if(q == 10){q10++;
              }else if(q == 21){q21++;
              }else if(q == 30){q30++;
              }else if(q == 32){q32++;
              }else if(q == 40){q40++;
              }else if(q == 41){q41++;
              }else if(q == 42){q42++;
              }else if(q == 43){q43++;
              }else if(q == 50){q50++;
              }else if(q == 51){q51++;
              }else if(q == 52){q52++;
              }else if(q == 53){q53++;
              }else if(q == 410){q410++;
              }else if(q == 421){q421++;
              }else if(q == 430){q430++;
              }else if(q == 432){q432++;
              }else if(q == 510){q510++;
              }else if(q == 521){q521++;
              }else if(q == 530){q530++;
              }else if(q == 532){q532++;
              }else{System.out.println("Error: Invlid Quadrant");}
           }
           //q10, q21, q30, q32, q40, q41, q42, q43, q50, q51, q52, q53, q410, q421, q430, q432, q510, q521, q530, q532;
           aveQ[0] = aveQ[0] + q0;   aveQ[1] = aveQ[1] + q1;   aveQ[2] = aveQ[2] + q2;
           aveQ[3] = aveQ[3] + q3;   aveQ[4] = aveQ[4] + q4;   aveQ[5] = aveQ[5] + q5;
           aveQ[6] = aveQ[6] + q10;   aveQ[7] = aveQ[7] + q21;   aveQ[8] = aveQ[8] + q30;
           aveQ[9] = aveQ[9] + q32;   aveQ[10] = aveQ[10] + q40;   aveQ[11] = aveQ[11] + q41;
           aveQ[12] = aveQ[12] + q42;   aveQ[13] = aveQ[13] + q43;   aveQ[14] = aveQ[14] + q50;
           aveQ[15] = aveQ[15] + q51;   aveQ[16] = aveQ[16] + q52;   aveQ[17] = aveQ[17] + q53;
           aveQ[18] = aveQ[18] + q410;   aveQ[19] = aveQ[19] + q421;   aveQ[20] = aveQ[20] + q430;
           aveQ[21] = aveQ[21] + q432;   aveQ[22] = aveQ[22] + q510;   aveQ[23] = aveQ[23] + q521;
           aveQ[24] = aveQ[24] + q530;   aveQ[25] = aveQ[25] + q532;
           total1 = q0 + q1 + q2 + q3 + q4 + q5 + q10 + q21 + q30 + q32 + q40 + q41 + q42 + q43 + q50 + q51 + q52 + q53 + q410 + q421 + q430 + q432 + q510 + q521 + q530 + q532;
           System.out.println("k "+k+" = "+q0+" "+q1+" "+q2+" "+q3+" "+q4+" "+q5+" "+q10+" "+q21+" "+q30+" "+q32+" "+q40+" "+q41+" "+q42+" "+q43+" "+q50+" "+q51+" "+q52+" "+q53+" "+q410+" "+q421+" "+q430+" "+q432+" "+q510+" "+q521+" "+q530+" "+q532+" = "+total1);
           q0 = q0 + (q10/2) + (q30/2) + (q40/2) + (q50/2) + (q410/3) + (q430/3) + (q510/3) + (q530/3);
           q1 = q1 + (q10/2) + (q21/2) + (q41/2) + (q51/2) + (q410/3) + (q421/3) + (q510/3) + (q521/3);
           q2 = q2 + (q21/2) + (q32/2) + (q42/2) + (q52/2) + (q421/3) + (q432/3) + (q521/3) + (q532/3);
           q3 = q3 + (q30/2) + (q32/2) + (q43/2) + (q53/2) + (q430/3) + (q432/3) + (q530/3) + (q532/3);
           q4 = q4 + (q40/2) + (q41/2) + (q42/2) + (q43/2) + (q410/3) + (q421/3) + (q430/3) + (q432/3);
           q5 = q5 + (q50/2) + (q51/2) + (q52/2) + (q53/2) + (q510/3) + (q521/3) + (q530/3) + (q532/3);           
           float aveDiff_SIX = 0;
           aveDiff_SIX = Math.abs((multi/6)-q0) + Math.abs((multi/6)-q1) + Math.abs((multi/6)-q2) + Math.abs((multi/6)-q3) + Math.abs((multi/6)-q4) + Math.abs((multi/6)-q5);
           aveDiff_TOTAL = aveDiff_TOTAL + (aveDiff_SIX/6);
           System.out.println("k "+k+" = "+q0+" "+q1+" "+q2+" "+q3+" "+q4+" "+q5+"    Average bin difference = "+(aveDiff_SIX/6)+" = "+((aveDiff_SIX/6)/(multi)*100)+"%");
           count++;
        }
        aveQ[0] /= count;  aveQ[1] /= count; aveQ[2] /= count; aveQ[3] /= count; aveQ[4] /= count; aveQ[5] /= count;
        aveQ[6] /= count;  aveQ[7] /= count; aveQ[8] /= count; aveQ[9] /= count; aveQ[10] /= count; aveQ[11] /= count;
        aveQ[12] /= count;  aveQ[13] /= count; aveQ[14] /= count; aveQ[15] /= count; aveQ[16] /= count; aveQ[17] /= count;
        aveQ[18] /= count;  aveQ[19] /= count; aveQ[20] /= count; aveQ[21] /= count; aveQ[22] /= count; aveQ[23] /= count;
        aveQ[24] /= count;  aveQ[25] /= count;
        total2 = aveQ[0] + aveQ[1] + aveQ[2] + aveQ[3] + aveQ[4] + aveQ[5] + aveQ[6] + aveQ[7] + aveQ[8] + aveQ[9] + aveQ[10] + aveQ[11] + aveQ[12];
        total2 = total2 + aveQ[13] + aveQ[14] + aveQ[15] + aveQ[16] + aveQ[17] + aveQ[18] + aveQ[19] + aveQ[20] + aveQ[21] + aveQ[22] + aveQ[23] + aveQ[24] + aveQ[25];
        System.out.print("ave = "+aveQ[0]+" "+aveQ[1]+" "+aveQ[2]+" "+aveQ[3]+" "+aveQ[4]+" "+aveQ[5]+" "+aveQ[6]+" "+aveQ[7]+" "+aveQ[8]+" "+aveQ[9]+" "+aveQ[10]+" "+aveQ[11]+" "+aveQ[12]+" ");
        System.out.println(aveQ[13]+" "+aveQ[14]+" "+aveQ[15]+" "+aveQ[16]+" "+aveQ[17]+" "+aveQ[18]+" "+aveQ[19]+" "+aveQ[20]+" "+aveQ[21]+" "+aveQ[22]+" "+aveQ[23]+" "+aveQ[24]+" "+aveQ[25]+" = "+total2);
        aveQ[0] = aveQ[0] + (aveQ[6]/2) + (aveQ[8]/2) + (aveQ[10]/2) + (aveQ[14]/2) + (aveQ[18]/3) + (aveQ[20]/3) + (aveQ[22]/3) + (aveQ[24]/3);    //ok
        aveQ[1] = aveQ[1] + (aveQ[6]/2) + (aveQ[7]/2) + (aveQ[11]/2) + (aveQ[15]/2) + (aveQ[18]/3) + (aveQ[19]/3) + (aveQ[22]/3) + (aveQ[23]/3);    //ok
        aveQ[2] = aveQ[2] + (aveQ[7]/2) + (aveQ[9]/2) + (aveQ[12]/2) + (aveQ[16]/2) + (aveQ[19]/3) + (aveQ[21]/3) + (aveQ[23]/3) + (aveQ[25]/3);    //ok
        aveQ[3] = aveQ[3] + (aveQ[8]/2) + (aveQ[9]/2) + (aveQ[13]/2) + (aveQ[17]/2) + (aveQ[20]/3) + (aveQ[21]/3) + (aveQ[24]/3) + (aveQ[25]/3);    //ok
        aveQ[4] = aveQ[4] + (aveQ[10]/2) + (aveQ[11]/2) + (aveQ[12]/2) + (aveQ[13]/2) + (aveQ[18]/3) + (aveQ[19]/3) + (aveQ[20]/3) + (aveQ[21]/3);  //ok
        aveQ[5] = aveQ[5] + (aveQ[14]/2) + (aveQ[15]/2) + (aveQ[16]/2) + (aveQ[17]/2) + (aveQ[22]/3) + (aveQ[23]/3) + (aveQ[24]/3) + (aveQ[25]/3);  //ok
        System.out.println("ave = "+aveQ[0]+" "+aveQ[1]+" "+aveQ[2]+" "+aveQ[3]+" "+aveQ[4]+" "+aveQ[5]+"  Total Average bin difference = "+(aveDiff_TOTAL/tests)+" = "+((aveDiff_TOTAL/tests)/(multi)*100)+"%");
        System.out.println();
        System.out.println("H = c0, c1, c2, c3, c4, c5, c6, c7");
        total1 = c0 + c1 + c2 + c3 + c4 + c5 + c6 + c7;
        System.out.println("C = "+c0+" "+c1+" "+c2+" "+c3+" "+c4+" "+c5+" "+c6+" "+c7+" = "+total1);
        System.out.println();
    }
    public void testWaterGeneratinQuadrantAllocation(TheMatrix TM){
        double rx, ry, rz;
        double[] rx2, ry2, rz2;
        TM.sitesMol = 4;
        rx2 = new double[TM.sitesMol];   ry2 = new double[TM.sitesMol];
        rz2 = new double[TM.sitesMol];
        double a1, a2, a3;
        double eAngx, eAngy, eAngz;
        int total1 = 0; int total2 = 0; int total3 = 0;
        double total1Av, total2Av, total3Av;
        total1Av = total2Av = total3Av = 0;
        int q, q0, q1, q2, q3, q4, q5, q10, q21, q30, q32, q40, q41, q42, q43, q50, q51, q52, q53, q410, q421, q430, q432, q510, q521, q530, q532;
        q = q0 = q1 = q2 = q3 = q4 = q5 = q10 = q21 = q30 = q32 = q40 = q41 = q42 = q43 = q50 = q51 = q52 = q53 = q410 = q421 = q430 = q432 = q510 = q521 = q530 = q532 = 0;
        int qq, qq0, qq1, qq2, qq3, qq4, qq5, qq10, qq21, qq30, qq32, qq40, qq41, qq42, qq43, qq50, qq51, qq52, qq53, qq410, qq421, qq430, qq432, qq510, qq521, qq530, qq532;
        qq = qq0 = qq1 = qq2 = qq3 = qq4 = qq5 = qq10 = qq21 = qq30 = qq32 = qq40 = qq41 = qq42 = qq43 = qq50 = qq51 = qq52 = qq53 = qq410 = qq421 = qq430 = qq432 = qq510 = qq521 = qq530 = qq532 = 0;
        int qqq, qqq0, qqq1, qqq2, qqq3, qqq4, qqq5, qqq10, qqq21, qqq30, qqq32, qqq40, qqq41, qqq42, qqq43, qqq50, qqq51, qqq52, qqq53, qqq410, qqq421, qqq430, qqq432, qqq510, qqq521, qqq530, qqq532;
        qqq = qqq0 = qqq1 = qqq2 = qqq3 = qqq4 = qqq5 = qqq10 = qqq21 = qqq30 = qqq32 = qqq40 = qqq41 = qqq42 = qqq43 = qqq50 = qqq51 = qqq52 = qqq53 = qqq410 = qqq421 = qqq430 = qqq432 = qqq510 = qqq521 = qqq530 = qqq532 = 0;
        int count = 0;
        double q_u1, q_u2, q_u3, q_u4;
        double[] tq = new double[4];
        double[] p = new double[10];
        int k, k1, k2;
        double s = 0;
        double tx, ty, tz, temp;
        double[] rMatT = new double[9];
        double[] rmx, rmy, rmz;
        rmx = new double[4]; rmy = new double[4]; rmz = new double[4];
        for(int j = 0; j < TM.sitesMol; j++){
            rmx[j] = 0.0; rmy[j] = 0.0; rmz[j] = 0.0;
        }       
        rmy[2] = 0.2402;   rmz[2] = 0.186;
        rmy[3] =-0.2402;   rmz[3] = 0.186;
        double[]  qAv = new double[26];
        double[] qqAv = new double[26];
        double[] qqqAv = new double[26];
        int multi = 256;
        float qaveDiff_TOTAL = 0;
        float qqaveDiff_TOTAL = 0;
        float qqqaveDiff_TOTAL = 0;
        qAv[0]  = qAv[1]  = qAv[2]  = qAv[3]  = qAv[4]  = qAv[5]  = qAv[6]  = qAv[7]  = qAv[8]  = qAv[9]  = qAv[10] = qAv[11] = qAv[12] = 0;
        qAv[13] = qAv[14] = qAv[15] = qAv[16] = qAv[17] = qAv[18] = qAv[19] = qAv[20] = qAv[21] = qAv[22] = qAv[23] = qAv[24] = qAv[25] = 0;
        qqAv[0]  = qqAv[1]  = qqAv[2]  = qqAv[3]  = qqAv[4]  = qqAv[5]  = qqAv[6]  = qqAv[7]  = qqAv[8]  = qqAv[9]  = qqAv[10] = qqAv[11] = qqAv[12] = 0;
        qqAv[13] = qqAv[14] = qqAv[15] = qqAv[16] = qqAv[17] = qqAv[18] = qqAv[19] = qqAv[20] = qqAv[21] = qqAv[22] = qqAv[23] = qqAv[24] = qqAv[25] = 0;
        qqqAv[0]  = qqqAv[1]  = qqqAv[2]  = qqqAv[3]  = qqqAv[4]  = qqqAv[5]  = qqqAv[6]  = qqqAv[7]  = qqqAv[8]  = qqqAv[9]  = qqqAv[10] = qqqAv[11] = qqqAv[12] = 0;
        qqqAv[13] = qqqAv[14] = qqqAv[15] = qqqAv[16] = qqqAv[17] = qqqAv[18] = qqqAv[19] = qqqAv[20] = qqqAv[21] = qqqAv[22] = qqqAv[23] = qqqAv[24] = qqqAv[25] = 0;
        TM.initMatrix();
        double dy, dz, distHHx, distHHy, distHHz;
        for(int kk = 0; kk < tests; kk++){
           q = q0 = q1 = q2 = q3 = q4 = q5 = q10 = q21 = q30 = q32 = q40 = q41 = q42 = q43 = q50 = q51 = q52 = q53 = q410 = q421 = q430 = q432 = q510 = q521 = q530 = q532 = 0;
           qq = qq0 = qq1 = qq2 = qq3 = qq4 = qq5 = qq10 = qq21 = qq30 = qq32 = qq40 = qq41 = qq42 = qq43 = qq50 = qq51 = qq52 = qq53 = qq410 = qq421 = qq430 = qq432 = qq510 = qq521 = qq530 = qq532 = 0;
           qqq = qqq0 = qqq1 = qqq2 = qqq3 = qqq4 = qqq5 = qqq10 = qqq21 = qqq30 = qqq32 = qqq40 = qqq41 = qqq42 = qqq43 = qqq50 = qqq51 = qqq52 = qqq53 = qqq410 = qqq421 = qqq430 = qqq432 = qqq510 = qqq521 = qqq530 = qqq532 = 0;
           for(int ll = 0; ll < multi; ll++){
              rx = Math.random();
              ry = Math.random();
              rz = Math.random();
              if(Math.random() < 0.5){rx=-1*rx;}
              if(Math.random() < 0.5){ry=-1*ry;}
              if(Math.random() < 0.5){rz=-1*rz;}
              q = getQuadrant(rx,ry,rz);
                 if(q == 0){q0++;
              }else if(q == 1){q1++;
              }else if(q == 2){q2++;
              }else if(q == 3){q3++;
              }else if(q == 4){q4++;
              }else if(q == 5){q5++;
              }else if(q == 10){q10++;
              }else if(q == 21){q21++;
              }else if(q == 30){q30++;
              }else if(q == 32){q32++;
              }else if(q == 40){q40++;
              }else if(q == 41){q41++;
              }else if(q == 42){q42++;
              }else if(q == 43){q43++;
              }else if(q == 50){q50++;
              }else if(q == 51){q51++;
              }else if(q == 52){q52++;
              }else if(q == 53){q53++;
              }else if(q == 410){q410++;
              }else if(q == 421){q421++;
              }else if(q == 430){q430++;
              }else if(q == 432){q432++;
              }else if(q == 510){q510++;
              }else if(q == 521){q521++;
              }else if(q == 530){q530++;
              }else if(q == 532){q532++;
              }else{System.out.println("Error: Invalid Quadrant");}
     //       System.out.println("(rx,ry,rz)=("+rx+","+ry+","+rz+")  q="+q);
              eAngx = Math.atan2(rx,ry);         //a1 = 2*eAngx*Math.PI;     //a1 = 0.5 * eAngy;
              eAngy = Math.acos(rz);             //a2 = 2*eAngy*Math.PI;     //a2 = 0.5 * (eAngx - eAngz);
              temp = Math.random();
         //   if(Math.random() < 0.5){temp=-1*temp;}
              eAngz = 2*Math.PI*temp;   //a3 = 2*eAngz*Math.PI;     //a3 = 0.5 * (eAngx + eAngz);
     //       System.out.println("(eAngx,eAngy,eAngz) = ("+eAngx+","+eAngy+","+eAngz+")");
              a1 = 0.5*eAngy;
              a2 = 0.5*(eAngx-eAngz);
              a3 = 0.5*(eAngx+eAngz);
   //           System.out.println("(a1,a2,a3) = ("+a1+","+a2+","+a3+")");
              q_u1 = Math.sin(a1) * Math.cos(a2);
              q_u2 = Math.sin(a1) * Math.sin(a2);
              q_u3 = Math.cos(a1) * Math.sin(a3);
              q_u4 = Math.cos(a1) * Math.cos(a3);
  //            System.out.println("(q1,q2,q3,q4) = ("+q_u1+","+q_u2+","+q_u3+","+q_u4+")");
              tq[0] = q_u1;  tq[1] = q_u2;  tq[2] = q_u3;    tq[3] = q_u4;
              for(k = 0, k2 = 0; k2 < 4; k2++){
                 for(k1 = k2; k1 < 4; k1++, k++){
                     p[k] = 2.0*tq[k1]*tq[k2];
                 }
              }
              rMatT[0] = p[0] + p[9] - 1;   rMatT[4] = p[4] + p[9] - 1;   rMatT[8] = p[7] + p[9] - 1;
              s = -1.0;    //Transpose = 1
              rMatT[1] = p[1] + s * p[8];   rMatT[3] = p[1] - s * p[8];   rMatT[2] = p[2] - s * p[6];
              rMatT[6] = p[2] + s * p[6];   rMatT[5] = p[5] + s * p[3];   rMatT[7] = p[5] - s * p[3];
              for(int j = 0; j < TM.sitesMol; j++){   // tip3p only
                 tx = rMatT[0]*rmx[j] + rMatT[3]*rmy[j] + rMatT[6]*rmz[j];
                 ty = rMatT[1]*rmx[j] + rMatT[4]*rmy[j] + rMatT[7]*rmz[j];
                 tz = rMatT[2]*rmx[j] + rMatT[5]*rmy[j] + rMatT[8]*rmz[j];
                 rx2[j] = tx;  ry2[j] = ty;   rz2[j] = tz;
              }              
              dy = (rmy[2]-rmy[3])/2;  dz = (rmz[2]-rmz[3])/2;
              distHHx = (rx2[2]+rx2[3])/2;    distHHy = (ry2[2]+ry2[3])/2;   distHHz = (rz2[2]+rz2[3])/2;
              /**
               * This section would just print some distances to make sure the water were in the right geometry of water.
               *
              if(ll < 10){                 
                 System.out.println(ll+" O("+rx2[0]+","+ry2[0]+","+rz2[0]+") H1("+rx2[2]+","+ry2[2]+","+rz2[2]+") H2("+rx2[3]+","+ry2[3]+","+rz2[3]+")");
                 System.out.println("     d(O,H1)="+(Math.sqrt(((rx2[0]-rx2[2])*(rx2[0]-rx2[2]))+((ry2[0]-ry2[2])*(ry2[0]-ry2[2]))+((rz2[0]-rz2[2])*(rz2[0]-rz2[2])))));
                 System.out.println("     d(O,H2)="+(Math.sqrt(((rx2[0]-rx2[3])*(rx2[0]-rx2[3]))+((ry2[0]-ry2[3])*(ry2[0]-ry2[3]))+((rz2[0]-rz2[3])*(rz2[0]-rz2[3])))));
                 System.out.println("    d(H1,H2)="+(Math.sqrt(((rx2[2]-rx2[3])*(rx2[2]-rx2[3]))+((ry2[2]-ry2[3])*(ry2[2]-ry2[3]))+((rz2[2]-rz2[3])*(rz2[2]-rz2[3])))));                 
                 System.out.println("    vect(x,y,z)=("+distHHx+","+distHHy+","+distHHz+")");
                 System.out.println("     d(V,O)="+(Math.sqrt(((distHHx-rx2[0])*(distHHx-rx2[0]))+((distHHy-ry2[0])*(distHHy-ry2[0]))+((distHHz-rz2[0])*(distHHz-rz2[0])))));
                 System.out.println("     d(V,H1)="+(Math.sqrt(((distHHx-rx2[2])*(distHHx-rx2[2]))+((distHHy-ry2[2])*(distHHy-ry2[2]))+((distHHz-rz2[2])*(distHHz-rz2[2])))));
                 System.out.println("     d(V,H2)="+(Math.sqrt(((distHHx-rx2[3])*(distHHx-rx2[3]))+((distHHy-ry2[3])*(distHHy-ry2[3]))+((distHHz-rz2[3])*(distHHz-rz2[3])))));
              }*/           
              total1 = q0 + q1 + q2 + q3 + q4 + q5 + q10 + q21 + q30 + q32 + q40 + q41 + q42 + q43 + q50 + q51 + q52 + q53 + q410 + q421 + q430 + q432 + q510 + q521 + q530 + q532;
              total1Av = total1Av + total1;
              qq = getQuadrant(distHHx,distHHy,distHHz);
                    if(qq == 0){qq0++;
              }else if(qq == 1){qq1++;
              }else if(qq == 2){qq2++;
              }else if(qq == 3){qq3++;
              }else if(qq == 4){qq4++;
              }else if(qq == 5){qq5++;
              }else if(qq == 10){qq10++;
              }else if(qq == 21){qq21++;
              }else if(qq == 30){qq30++;
              }else if(qq == 32){qq32++;
              }else if(qq == 40){qq40++;
              }else if(qq == 41){qq41++;
              }else if(qq == 42){qq42++;
              }else if(qq == 43){qq43++;
              }else if(qq == 50){qq50++;
              }else if(qq == 51){qq51++;
              }else if(qq == 52){qq52++;
              }else if(qq == 53){qq53++;
              }else if(qq == 410){qq410++;
              }else if(qq == 421){qq421++;
              }else if(qq == 430){qq430++;
              }else if(qq == 432){qq432++;
              }else if(qq == 510){qq510++;
              }else if(qq == 521){qq521++;
              }else if(qq == 530){qq530++;
              }else if(qq == 532){qq532++;
              }else{System.out.println("Error: Invalid Quadrant");}
              total2 = qq0 + qq1 + qq2 + qq3 + qq4 + qq5 + qq10 + qq21 + qq30 + qq32 + qq40 + qq41 + qq42 + qq43 + qq50 + qq51 + qq52 + qq53 + qq410 + qq421 + qq430 + qq432 + qq510 + qq521 + qq530 + qq532;
              total2Av = total2Av + total2;
              double[] mc, mt;
              //    double tx, ty, tz;
              
              mc = new double[9];          mt = new double[9];
              TM.wvx[0] = 0.0;  TM.wvy[0] = 0.0; TM.wvz[0] = 0.1;
              for(int n = 0; n < 50; n++){
                  tx = Math.random()*Math.PI*2;
                  ty = Math.random()*Math.PI*2;
                  tz = Math.random()*Math.PI*2;
                  BuildStepRmat(mc, tx, ty, tz);
                  MulMat(mt, mc, rMatT, 3, 0);
                  rMatT[0] = mt[0];    rMatT[1] = mt[1];    rMatT[2] = mt[2];
                  rMatT[3] = mt[3];    rMatT[4] = mt[4];    rMatT[5] = mt[5];
                  rMatT[6] = mt[6];    rMatT[7] = mt[7];    rMatT[8] = mt[8];
                  for(int j = 0; j < TM.sitesMol; j++){   // tip3p only
                      tx = rMatT[0]*rmx[j] + rMatT[3]*rmy[j] + rMatT[6]*rmz[j];
                      ty = rMatT[1]*rmx[j] + rMatT[4]*rmy[j] + rMatT[7]*rmz[j];
                      tz = rMatT[2]*rmx[j] + rMatT[5]*rmy[j] + rMatT[8]*rmz[j];
                      rx2[j] = tx;  ry2[j] = ty;   rz2[j] = tz;
                 }
              }
              dy = (rmy[2]-rmy[3])/2;  dz = (rmz[2]-rmz[3])/2;
              distHHx = (rx2[2]+rx2[3])/2;    distHHy = (ry2[2]+ry2[3])/2;   distHHz = (rz2[2]+rz2[3])/2;
              qqq = getQuadrant(distHHx,distHHy,distHHz);
                    if(qqq == 0){qqq0++;
              }else if(qqq == 1){qqq1++;
              }else if(qqq == 2){qqq2++;
              }else if(qqq == 3){qqq3++;
              }else if(qqq == 4){qqq4++;
              }else if(qqq == 5){qqq5++;
              }else if(qqq == 10){qqq10++;
              }else if(qqq == 21){qqq21++;
              }else if(qqq == 30){qqq30++;
              }else if(qqq == 32){qqq32++;
              }else if(qqq == 40){qqq40++;
              }else if(qqq == 41){qqq41++;
              }else if(qqq == 42){qqq42++;
              }else if(qqq == 43){qqq43++;
              }else if(qqq == 50){qqq50++;
              }else if(qqq == 51){qqq51++;
              }else if(qqq == 52){qqq52++;
              }else if(qqq == 53){qqq53++;
              }else if(qqq == 410){qqq410++;
              }else if(qqq == 421){qqq421++;
              }else if(qqq == 430){qqq430++;
              }else if(qqq == 432){qqq432++;
              }else if(qqq == 510){qqq510++;
              }else if(qqq == 521){qqq521++;
              }else if(qqq == 530){qqq530++;
              }else if(qqq == 532){qqq532++;
              }else{System.out.println("Error: Invalid Quadrant");}
              total3 = qqq0 + qqq1 + qqq2 + qqq3 + qqq4 + qqq5 + qqq10 + qqq21 + qqq30 + qqq32 + qqq40 + qqq41 + qqq42 + qqq43 + qqq50 + qqq51 + qqq52 + qqq53 + qqq410 + qqq421 + qqq430 + qqq432 + qqq510 + qqq521 + qqq530 + qqq532;
              total3Av = total3Av + total3;
           }
           qAv[0]  = qAv[0]  + q0;   qAv[1]  = qAv[1]  + q1;   qAv[2]  = qAv[2]  + q2;   qAv[3]  = qAv[3]  + q3;   qAv[4]  = qAv[4]  + q4;   qAv[5]  = qAv[5]  + q5;
           qAv[6]  = qAv[6]  + q10;  qAv[7]  = qAv[7]  + q21;  qAv[8]  = qAv[8]  + q30;  qAv[9]  = qAv[9]  + q32;  qAv[10] = qAv[10] + q40;  qAv[11] = qAv[11] + q41;
           qAv[12] = qAv[12] + q42;  qAv[13] = qAv[13] + q43;  qAv[14] = qAv[14] + q50;  qAv[15] = qAv[15] + q51;  qAv[16] = qAv[16] + q52;  qAv[17] = qAv[17] + q53;
           qAv[18] = qAv[18] + q410; qAv[19] = qAv[19] + q421; qAv[20] = qAv[20] + q430; qAv[21] = qAv[21] + q432; qAv[22] = qAv[22] + q510; qAv[23] = qAv[23] + q521;
           qAv[24] = qAv[24] + q530; qAv[25] = qAv[25] + q532;
      //     System.out.println(kk+" V "+q0+" + "+q1+" + "+q2+" + "+q3+" + "+q4+" + "+q5+" + "+q10+" + "+q21+" + "+q30+" + "+q32+" + "+q40+" + "+q41+" + "+q42+" + "+q43+" + "+q50+" + "+q51+" + "+q52+" + "+q53+" + "+q410+" + "+q421+" + "+q430+" + "+q432+" + "+q510+" + "+q521+" + "+q530+" + "+q532+" = "+total1);
           q0 = q0 + (q10/2) + (q30/2) + (q40/2) + (q50/2) + (q410/3) + (q430/3) + (q510/3) + (q530/3);
           q1 = q1 + (q10/2) + (q21/2) + (q41/2) + (q51/2) + (q410/3) + (q421/3) + (q510/3) + (q521/3);
           q2 = q2 + (q21/2) + (q32/2) + (q42/2) + (q52/2) + (q421/3) + (q432/3) + (q521/3) + (q532/3);
           q3 = q3 + (q30/2) + (q32/2) + (q43/2) + (q53/2) + (q430/3) + (q432/3) + (q530/3) + (q532/3);
           q4 = q4 + (q40/2) + (q41/2) + (q42/2) + (q43/2) + (q410/3) + (q421/3) + (q430/3) + (q432/3);
           q5 = q5 + (q50/2) + (q51/2) + (q52/2) + (q53/2) + (q510/3) + (q521/3) + (q530/3) + (q532/3);
           float qaveDiff_SIX = 0;
           qaveDiff_SIX = Math.abs((multi/6)-q0) + Math.abs((multi/6)-q1) + Math.abs((multi/6)-q2) + Math.abs((multi/6)-q3) + Math.abs((multi/6)-q4) + Math.abs((multi/6)-q5);
           qaveDiff_TOTAL = qaveDiff_TOTAL + (qaveDiff_SIX/6);
           System.out.println(kk+" V "+q0+" "+q1+" "+q2+" "+q3+" "+q4+" "+q5+"    Average bin difference = "+(qaveDiff_SIX/6)+" = "+((qaveDiff_SIX/6)/(multi)*100)+"%");
           qqAv[0]  = qqAv[0]  + qq0;   qqAv[1]  = qqAv[1]  + qq1;   qqAv[2]  = qqAv[2]  + qq2;   qqAv[3]  = qqAv[3]  + qq3;   qqAv[4]  = qqAv[4]  + qq4;   qqAv[5]  = qqAv[5] +  qq5;
           qqAv[6]  = qqAv[6]  + qq10;  qqAv[7]  = qqAv[7]  + qq21;  qqAv[8]  = qqAv[8]  + qq30;  qqAv[9]  = qqAv[9]  + qq32;  qqAv[10] = qqAv[10] + qq40;  qqAv[11] = qqAv[11] + qq41;
           qqAv[12] = qqAv[12] + qq42;  qqAv[13] = qqAv[13] + qq43;  qqAv[14] = qqAv[14] + qq50;  qqAv[15] = qqAv[15] + qq51;  qqAv[16] = qqAv[16] + qq52;  qqAv[17] = qqAv[17] + qq53;
           qqAv[18] = qqAv[18] + qq410; qqAv[19] = qqAv[19] + qq421; qqAv[20] = qqAv[20] + qq430; qqAv[21] = qqAv[21] + qq432; qqAv[22] = qqAv[22] + qq510; qqAv[23] = qqAv[23] + qq521;
           qqAv[24] = qqAv[24] + qq530; qqAv[25] = qqAv[25] + qq532;           
 //          System.out.println(kk+" W "+qq0+" + "+qq1+" + "+qq2+" + "+qq3+" + "+qq4+" + "+qq5+" + "+qq10+" + "+qq21+" + "+qq30+" + "+qq32+" + "+qq40+" + "+qq41+" + "+qq42+" + "+qq43+" + "+qq50+" + "+qq51+" + "+qq52+" + "+qq53+" + "+qq410+" + "+qq421+" + "+qq430+" + "+qq432+" + "+qq510+" + "+qq521+" + "+qq530+" + "+qq532+" = "+total2);
           qq0 = qq0 + (qq10/2) + (qq30/2) + (qq40/2) + (qq50/2) + (qq410/3) + (qq430/3) + (qq510/3) + (qq530/3);
           qq1 = qq1 + (qq10/2) + (qq21/2) + (qq41/2) + (qq51/2) + (qq410/3) + (qq421/3) + (qq510/3) + (qq521/3);
           qq2 = qq2 + (qq21/2) + (qq32/2) + (qq42/2) + (qq52/2) + (qq421/3) + (qq432/3) + (qq521/3) + (qq532/3);
           qq3 = qq3 + (qq30/2) + (qq32/2) + (qq43/2) + (qq53/2) + (qq430/3) + (qq432/3) + (qq530/3) + (qq532/3);
           qq4 = qq4 + (qq40/2) + (qq41/2) + (qq42/2) + (qq43/2) + (qq410/3) + (qq421/3) + (qq430/3) + (qq432/3);
           qq5 = qq5 + (qq50/2) + (qq51/2) + (qq52/2) + (qq53/2) + (qq510/3) + (qq521/3) + (qq530/3) + (qq532/3);
           float qqaveDiff_SIX = 0;
           qqaveDiff_SIX = Math.abs((multi/6)-qq0) + Math.abs((multi/6)-qq1) + Math.abs((multi/6)-qq2) + Math.abs((multi/6)-qq3) + Math.abs((multi/6)-qq4) + Math.abs((multi/6)-qq5);
           qqaveDiff_TOTAL = qqaveDiff_TOTAL + (qqaveDiff_SIX/6);
           System.out.println(kk+" W "+qq0+" "+qq1+" "+qq2+" "+qq3+" "+qq4+" "+qq5+"    Average bin difference = "+(qqaveDiff_SIX/6)+" = "+((qqaveDiff_SIX/6)/(multi)*100)+"%");
           
           qqqAv[0]  = qqqAv[0]  + qqq0;   qqqAv[1]  = qqqAv[1]  + qqq1;   qqqAv[2]  = qqqAv[2]  + qqq2;   qqqAv[3]  = qqqAv[3]  + qqq3;   qqqAv[4]  = qqqAv[4]  + qqq4;   qqqAv[5]  = qqqAv[5] +  qqq5;
           qqqAv[6]  = qqqAv[6]  + qqq10;  qqqAv[7]  = qqqAv[7]  + qqq21;  qqqAv[8]  = qqqAv[8]  + qqq30;  qqqAv[9]  = qqqAv[9]  + qqq32;  qqqAv[10] = qqqAv[10] + qqq40;  qqqAv[11] = qqqAv[11] + qqq41;
           qqqAv[12] = qqqAv[12] + qqq42;  qqqAv[13] = qqqAv[13] + qqq43;  qqqAv[14] = qqqAv[14] + qqq50;  qqqAv[15] = qqqAv[15] + qqq51;  qqqAv[16] = qqqAv[16] + qqq52;  qqqAv[17] = qqqAv[17] + qqq53;
           qqqAv[18] = qqqAv[18] + qqq410; qqqAv[19] = qqqAv[19] + qqq421; qqqAv[20] = qqqAv[20] + qqq430; qqqAv[21] = qqqAv[21] + qqq432; qqqAv[22] = qqqAv[22] + qqq510; qqqAv[23] = qqqAv[23] + qqq521;
           qqqAv[24] = qqqAv[24] + qqq530; qqqAv[25] = qqqAv[25] + qqq532;           
//           System.out.println(kk+" X "+qqq0+" + "+qqq1+" + "+qqq2+" + "+qqq3+" + "+qqq4+" + "+qqq5+" + "+qqq10+" + "+qqq21+" + "+qqq30+" + "+qqq32+" + "+qqq40+" + "+qqq41+" + "+qqq42+" + "+qqq43+" + "+qqq50+" + "+qqq51+" + "+qqq52+" + "+qqq53+" + "+qqq410+" + "+qqq421+" + "+qqq430+" + "+qqq432+" + "+qqq510+" + "+qqq521+" + "+qqq530+" + "+qqq532+" = "+total3);
           qqq0 = qqq0 + (qqq10/2) + (qqq30/2) + (qqq40/2) + (qqq50/2) + (qqq410/3) + (qqq430/3) + (qqq510/3) + (qqq530/3);
           qqq1 = qqq1 + (qqq10/2) + (qqq21/2) + (qqq41/2) + (qqq51/2) + (qqq410/3) + (qqq421/3) + (qqq510/3) + (qqq521/3);
           qqq2 = qqq2 + (qqq21/2) + (qqq32/2) + (qqq42/2) + (qqq52/2) + (qqq421/3) + (qqq432/3) + (qqq521/3) + (qqq532/3);
           qqq3 = qqq3 + (qqq30/2) + (qqq32/2) + (qqq43/2) + (qqq53/2) + (qqq430/3) + (qqq432/3) + (qqq530/3) + (qqq532/3);
           qqq4 = qqq4 + (qqq40/2) + (qqq41/2) + (qqq42/2) + (qqq43/2) + (qqq410/3) + (qqq421/3) + (qqq430/3) + (qqq432/3);
           qqq5 = qqq5 + (qqq50/2) + (qqq51/2) + (qqq52/2) + (qqq53/2) + (qqq510/3) + (qqq521/3) + (qqq530/3) + (qqq532/3);
           float qqqaveDiff_SIX = 0;
           qqqaveDiff_SIX = Math.abs((multi/6)-qqq0) + Math.abs((multi/6)-qqq1) + Math.abs((multi/6)-qqq2) + Math.abs((multi/6)-qqq3) + Math.abs((multi/6)-qqq4) + Math.abs((multi/6)-qqq5);
           qqqaveDiff_TOTAL = qqqaveDiff_TOTAL + (qqqaveDiff_SIX/6);
           System.out.println(kk+" X "+qqq0+" "+qqq1+" "+qqq2+" "+qqq3+" "+qqq4+" "+qqq5+"    Average bin difference = "+(qqqaveDiff_SIX/6)+" = "+((qqqaveDiff_SIX/6)/(multi)*100)+"%");
           count++;
        }
        qAv[0]  /= count;  qAv[1]  /= count; qAv[2]  /= count; qAv[3]  /= count; qAv[4]  /= count; qAv[5]  /= count;
        qAv[6]  /= count;  qAv[7]  /= count; qAv[8]  /= count; qAv[9]  /= count; qAv[10] /= count; qAv[11] /= count;
        qAv[12] /= count;  qAv[13] /= count; qAv[14] /= count; qAv[15] /= count; qAv[16] /= count; qAv[17] /= count;
        qAv[18] /= count;  qAv[19] /= count; qAv[20] /= count; qAv[21] /= count; qAv[22] /= count; qAv[23] /= count;
        qAv[24] /= count;  qAv[25] /= count;        
        total1Av = qAv[0]+qAv[1]+qAv[2]+qAv[3]+qAv[4]+qAv[5]+qAv[6]+qAv[7]+qAv[8]+qAv[9]+qAv[10]+qAv[11]+qAv[12];
        total1Av = total1Av+qAv[13]+qAv[14]+qAv[15]+qAv[16]+qAv[17]+qAv[18]+qAv[19]+qAv[20]+qAv[21]+qAv[22]+qAv[23]+qAv[24]+qAv[25];
        System.out.print("ave V = "+qAv[0]+" "+qAv[1]+" "+qAv[2]+" "+qAv[3]+" "+qAv[4]+" "+qAv[5]+" "+qAv[6]+" "+qAv[7]+" "+qAv[8]+" "+qAv[9]+" "+qAv[10]+" "+qAv[11]+" "+qAv[12]+" ");
        System.out.println(qAv[13]+" "+qAv[14]+" "+qAv[15]+" "+qAv[16]+" "+qAv[17]+" "+qAv[18]+" "+qAv[19]+" "+qAv[20]+" "+qAv[21]+" "+qAv[22]+" "+qAv[23]+" "+qAv[24]+" "+qAv[25]+" = "+total1Av);
        qAv[0] = qAv[0] + (qAv[6]/2) + (qAv[8]/2) + (qAv[10]/2) + (qAv[14]/2) + (qAv[18]/3) + (qAv[20]/3) + (qAv[22]/3) + (qAv[24]/3);    //ok
        qAv[1] = qAv[1] + (qAv[6]/2) + (qAv[7]/2) + (qAv[11]/2) + (qAv[15]/2) + (qAv[18]/3) + (qAv[19]/3) + (qAv[22]/3) + (qAv[23]/3);    //ok
        qAv[2] = qAv[2] + (qAv[7]/2) + (qAv[9]/2) + (qAv[12]/2) + (qAv[16]/2) + (qAv[19]/3) + (qAv[21]/3) + (qAv[23]/3) + (qAv[25]/3);    //ok
        qAv[3] = qAv[3] + (qAv[8]/2) + (qAv[9]/2) + (qAv[13]/2) + (qAv[17]/2) + (qAv[20]/3) + (qAv[21]/3) + (qAv[24]/3) + (qAv[25]/3);    //ok
        qAv[4] = qAv[4] + (qAv[10]/2) + (qAv[11]/2) + (qAv[12]/2) + (qAv[13]/2) + (qAv[18]/3) + (qAv[19]/3) + (qAv[20]/3) + (qAv[21]/3);  //ok
        qAv[5] = qAv[5] + (qAv[14]/2) + (qAv[15]/2) + (qAv[16]/2) + (qAv[17]/2) + (qAv[22]/3) + (qAv[23]/3) + (qAv[24]/3) + (qAv[25]/3);  //ok
        System.out.println("ave V = "+qAv[0]+" "+qAv[1]+" "+qAv[2]+" "+qAv[3]+" "+qAv[4]+" "+qAv[5]+"  Total Average bin difference = "+(qaveDiff_TOTAL/tests)+" = "+((qaveDiff_TOTAL/tests)/(multi)*100)+"%");
        qqAv[0]  /= count;  qqAv[1]  /= count; qqAv[2]  /= count; qqAv[3]  /= count; qqAv[4]  /= count; qqAv[5]  /= count;
        qqAv[6]  /= count;  qqAv[7]  /= count; qqAv[8]  /= count; qqAv[9]  /= count; qqAv[10] /= count; qqAv[11] /= count;
        qqAv[12] /= count;  qqAv[13] /= count; qqAv[14] /= count; qqAv[15] /= count; qqAv[16] /= count; qqAv[17] /= count;
        qqAv[18] /= count;  qqAv[19] /= count; qqAv[20] /= count; qqAv[21] /= count; qqAv[22] /= count; qqAv[23] /= count;
        qqAv[24] /= count;  qqAv[25] /= count;
        total2Av = qqAv[0]+qqAv[1]+qqAv[2]+qqAv[3]+qqAv[4]+qqAv[5]+qqAv[6]+qqAv[7]+qqAv[8]+qqAv[9]+qqAv[10]+qqAv[11]+qqAv[12];
        total2Av = total2Av+qqAv[13]+qqAv[14]+qqAv[15]+qqAv[16]+qqAv[17]+qqAv[18]+qqAv[19]+qqAv[20]+qqAv[21]+qqAv[22]+qqAv[23]+qqAv[24]+qqAv[25];
        System.out.print("ave W = "+qqAv[0]+" "+qqAv[1]+" "+qqAv[2]+" "+qqAv[3]+" "+qqAv[4]+" "+qqAv[5]+" "+qqAv[6]+" "+qqAv[7]+" "+qqAv[8]+" "+qqAv[9]+" "+qqAv[10]+" "+qqAv[11]+" "+qqAv[12]+" ");
        System.out.println(qqAv[13]+" "+qqAv[14]+" "+qqAv[15]+" "+qqAv[16]+" "+qqAv[17]+" "+qqAv[18]+" "+qqAv[19]+" "+qqAv[20]+" "+qqAv[21]+" "+qqAv[22]+" "+qqAv[23]+" "+qqAv[24]+" "+qqAv[25]+" = "+total2Av);
        qqAv[0] = qqAv[0] + (qqAv[6]/2) + (qqAv[8]/2) + (qqAv[10]/2) + (qqAv[14]/2) + (qqAv[18]/3) + (qqAv[20]/3) + (qqAv[22]/3) + (qqAv[24]/3);    //ok
        qqAv[1] = qqAv[1] + (qqAv[6]/2) + (qqAv[7]/2) + (qqAv[11]/2) + (qqAv[15]/2) + (qqAv[18]/3) + (qqAv[19]/3) + (qqAv[22]/3) + (qqAv[23]/3);    //ok
        qqAv[2] = qqAv[2] + (qqAv[7]/2) + (qqAv[9]/2) + (qqAv[12]/2) + (qqAv[16]/2) + (qqAv[19]/3) + (qqAv[21]/3) + (qqAv[23]/3) + (qqAv[25]/3);    //ok
        qqAv[3] = qqAv[3] + (qqAv[8]/2) + (qqAv[9]/2) + (qqAv[13]/2) + (qqAv[17]/2) + (qqAv[20]/3) + (qqAv[21]/3) + (qqAv[24]/3) + (qqAv[25]/3);    //ok
        qqAv[4] = qqAv[4] + (qqAv[10]/2) + (qqAv[11]/2) + (qqAv[12]/2) + (qqAv[13]/2) + (qqAv[18]/3) + (qqAv[19]/3) + (qqAv[20]/3) + (qqAv[21]/3);  //ok
        qqAv[5] = qqAv[5] + (qqAv[14]/2) + (qqAv[15]/2) + (qqAv[16]/2) + (qqAv[17]/2) + (qqAv[22]/3) + (qqAv[23]/3) + (qqAv[24]/3) + (qqAv[25]/3);  //ok
        System.out.println("ave W = "+qqAv[0]+" "+qqAv[1]+" "+qqAv[2]+" "+qqAv[3]+" "+qqAv[4]+" "+qqAv[5]+"  Total Average bin difference = "+(qqaveDiff_TOTAL/tests)+" = "+((qqaveDiff_TOTAL/tests)/(multi)*100)+"%");
        
        qqqAv[0]  /= count;  qqqAv[1]  /= count; qqqAv[2]  /= count; qqqAv[3]  /= count; qqqAv[4]  /= count; qqqAv[5]  /= count;
        qqqAv[6]  /= count;  qqqAv[7]  /= count; qqqAv[8]  /= count; qqqAv[9]  /= count; qqqAv[10] /= count; qqqAv[11] /= count;
        qqqAv[12] /= count;  qqqAv[13] /= count; qqqAv[14] /= count; qqqAv[15] /= count; qqqAv[16] /= count; qqqAv[17] /= count;
        qqqAv[18] /= count;  qqqAv[19] /= count; qqqAv[20] /= count; qqqAv[21] /= count; qqqAv[22] /= count; qqqAv[23] /= count;
        qqqAv[24] /= count;  qqqAv[25] /= count;
        total3Av = qqqAv[0]+qqqAv[1]+qqqAv[2]+qqqAv[3]+qqqAv[4]+qqqAv[5]+qqqAv[6]+qqqAv[7]+qqqAv[8]+qqqAv[9]+qqqAv[10]+qqqAv[11]+qqqAv[12];
        total3Av = total3Av+qqqAv[13]+qqqAv[14]+qqqAv[15]+qqqAv[16]+qqqAv[17]+qqqAv[18]+qqqAv[19]+qqqAv[20]+qqqAv[21]+qqqAv[22]+qqqAv[23]+qqqAv[24]+qqqAv[25];
        System.out.print("ave X = "+qqqAv[0]+" "+qqqAv[1]+" "+qqqAv[2]+" "+qqqAv[3]+" "+qqqAv[4]+" "+qqqAv[5]+" "+qqqAv[6]+" "+qqqAv[7]+" "+qqqAv[8]+" "+qqqAv[9]+" "+qqqAv[10]+" "+qqqAv[11]+" "+qqqAv[12]+" ");
        System.out.println(qqqAv[13]+" "+qqqAv[14]+" "+qqqAv[15]+" "+qqqAv[16]+" "+qqqAv[17]+" "+qqqAv[18]+" "+qqqAv[19]+" "+qqqAv[20]+" "+qqqAv[21]+" "+qqqAv[22]+" "+qqqAv[23]+" "+qqqAv[24]+" "+qqqAv[25]+" = "+total3Av);
        qqqAv[0] = qqqAv[0] + (qqqAv[6]/2) + (qqqAv[8]/2) + (qqqAv[10]/2) + (qqqAv[14]/2) + (qqqAv[18]/3) + (qqqAv[20]/3) + (qqqAv[22]/3) + (qqqAv[24]/3);    //ok
        qqqAv[1] = qqqAv[1] + (qqqAv[6]/2) + (qqqAv[7]/2) + (qqqAv[11]/2) + (qqqAv[15]/2) + (qqqAv[18]/3) + (qqqAv[19]/3) + (qqqAv[22]/3) + (qqqAv[23]/3);    //ok
        qqqAv[2] = qqqAv[2] + (qqqAv[7]/2) + (qqqAv[9]/2) + (qqqAv[12]/2) + (qqqAv[16]/2) + (qqqAv[19]/3) + (qqqAv[21]/3) + (qqqAv[23]/3) + (qqqAv[25]/3);    //ok
        qqqAv[3] = qqqAv[3] + (qqqAv[8]/2) + (qqqAv[9]/2) + (qqqAv[13]/2) + (qqqAv[17]/2) + (qqqAv[20]/3) + (qqqAv[21]/3) + (qqqAv[24]/3) + (qqqAv[25]/3);    //ok
        qqqAv[4] = qqqAv[4] + (qqqAv[10]/2) + (qqqAv[11]/2) + (qqqAv[12]/2) + (qqqAv[13]/2) + (qqqAv[18]/3) + (qqqAv[19]/3) + (qqqAv[20]/3) + (qqqAv[21]/3);  //ok
        qqqAv[5] = qqqAv[5] + (qqqAv[14]/2) + (qqqAv[15]/2) + (qqqAv[16]/2) + (qqqAv[17]/2) + (qqqAv[22]/3) + (qqqAv[23]/3) + (qqqAv[24]/3) + (qqqAv[25]/3);  //ok
        System.out.println("ave X = "+qqqAv[0]+" "+qqqAv[1]+" "+qqqAv[2]+" "+qqqAv[3]+" "+qqqAv[4]+" "+qqqAv[5]+"  Total Average bin difference = "+(qqqaveDiff_TOTAL/tests)+" = "+((qqqaveDiff_TOTAL/tests)/(multi)*100)+"%");
        System.out.println();
    }
    public float Round(float Rval, int Rpl){
           float p = (float)Math.pow(10,Rpl);
           Rval = Rval * p;
           float tmp = Math.round(Rval);
           return (float)tmp/p;
    }
    public int getQuadrant(double x, double y, double z){
        double ax, ay, az;
        ax = Math.abs(Round((float)x,3));   ay = Math.abs(Round((float)y,3));   az = Math.abs(Round((float)z,3));
              if((x>0)&&(ax>ay)&&(ax>az)){return 0;           // Volumen
        }else if((y>0)&&(ay>ax)&&(ay>az)){return 1;           // Volumen
        }else if((x<0)&&(ax>ay)&&(ax>az)){return 2;           // Volumen
        }else if((y<0)&&(ay>ax)&&(ay>az)){return 3;           // Volumen
        }else if((z>0)&&(az>ax)&&(az>ay)){return 4;           // Volumen
        }else if((z<0)&&(az>ax)&&(az>ay)){return 5;           // Volumen
        }else if((x>0)&&(y>0)&&(ax==ay)&&(ax>az)){return 10;  // Surface
        }else if((x>0)&&(y>0)&&(ax==ay)&&(ax<az)){
              if(z>0){return 4;}else if(z<0){return 5;}else{return -1;}       // Surface
        }else if((x<0)&&(y>0)&&(ax==ay)&&(ax>az)){return 21;  // Surface
        }else if((x<0)&&(y>0)&&(ax==ay)&&(ax<az)){
              if(z>0){return 4;}else if(z<0){return 5;}else{return -1;}       // Surface
        }else if((x<0)&&(y<0)&&(ax==ay)&&(ax>az)){return 32;  // Surface
        }else if((x<0)&&(y<0)&&(ax==ay)&&(ax<az)){
              if(z>0){return 4;}else if(z<0){return 5;}else{return -1;}       // Surface
        }else if((x>0)&&(y<0)&&(ax==ay)&&(ax>az)){return 30;  // Surface
        }else if((x>0)&&(y<0)&&(ax==ay)&&(ax<az)){
              if(z>0){return 4;}else if(x<0){return 5;}else{return -1;}       // Surface
        }else if((z>0)&&(x>0)&&(az==ax)&&(az>ay)){return 40;   // Surface
        }else if((z>0)&&(x>0)&&(az==ax)&&(az<ay)){
              if(y>0){return 1;}else if(y<0){return 3;}else{return -1;}       // Surface
        }else if((z>0)&&(y>0)&&(az==ay)&&(az>ax)){return 41;   // Surface
        }else if((z>0)&&(y>0)&&(az==ay)&&(az<ax)){
              if(x>0){return 0;}else if(x<0){return 2;}else{return -1;}       // Surface
        }else if((z>0)&&(x<0)&&(az==ax)&&(az>ay)){return 42;   // Surface
        }else if((z>0)&&(x<0)&&(az==ax)&&(az<ay)){
              if(y>0){return 1;}else if(y<0){return 3;}else{return -1;}       // Surface
        }else if((z>0)&&(y<0)&&(az==ay)&&(az>ax)){return 43;   // Surface
        }else if((z>0)&&(y<0)&&(az==ay)&&(az<ax)){
              if(x>0){return 0;}else if(x<0){return 2;}else{return -1;}       // Surface
        }else if((z<0)&&(x>0)&&(az==ax)&&(az>ay)){return 50;   // Surface
        }else if((z<0)&&(x>0)&&(az==ax)&&(az<ay)){
              if(y>0){return 1;}else if(y<0){return 3;}else{return -1;}       // Surface
        }else if((z<0)&&(y>0)&&(az==ay)&&(az>ax)){return 51;   // Surface
        }else if((z<0)&&(y>0)&&(az==ay)&&(az<ax)){
              if(x>0){return 0;}else if(x<0){return 2;}else{return -1;}       // Surface
        }else if((z<0)&&(x<0)&&(az==ax)&&(az>ay)){return 52;   // Surface
        }else if((z<0)&&(x<0)&&(az==ax)&&(az<ay)){
              if(y>0){return 1;}else if(y<0){return 3;}else{return -1;}       // Surface
        }else if((z<0)&&(y<0)&&(az==ay)&&(az>ax)){return 53;   // Surface
        }else if((z<0)&&(y<0)&&(az==ay)&&(az<ax)){
              if(x>0){return 0;}else if(x<0){return 2;}else{return -1;}      // Surface
        }else if((ax==ay)&&(ax==az)&&(ay==az)){
              if((x>0)&&(y>0)&&(z>0)){ return 410;
              }else if((x>0)&&(y>0)&&(z<0)){return 510;
              }else if((x>0)&&(y<0)&&(z>0)){return 430;
              }else if((x>0)&&(y<0)&&(z<0)){return 530;
              }else if((x<0)&&(y>0)&&(z>0)){return 421;
              }else if((x<0)&&(y>0)&&(z<0)){return 521;
              }else if((x<0)&&(y<0)&&(z>0)){return 432;
              }else if((x<0)&&(y<0)&&(z<0)){return 432;
              }else{return -1;}
        }else{
            return -1;
        }
    }
    public int getQuadrant1(double x, double y, double z){
        int Q = -1;
        double ax, ay, az;
        counter++;
        ax = Math.abs(Round((float)x,3));   ay = Math.abs(Round((float)y,3));   az = Math.abs(Round((float)z,3));
                     if(x >= 0 && y >= 0 && z >= 0){  c0++;  // 0
                              if(ax > ay && ax > az && ay > az){return 0;    // 0  (3,2,1)
                        }else if(ax > ay && ax > az && ay < az){return 0;    // 1  (3,1,2)
                        }else if(ax > ay && ax < az && ay < az){return 4;    // 3  (,,,)
                        }else if(ax < ay && ax > az && ay > az){return 1;    // 4  (,,,)
                        }else if(ax < ay && ax < az && ay > az){return 1;    // 6  (,,,)
                        }else if(ax < ay && ax < az && ay < az){return 4;    // 7  (,,,)
                        }else{
                                    if(ax == ay && ax == az && ay == az){return 0;
                          //    }else if(ax == ay && ax == az && ay != az){                                     
                          //    }else if(ax == ay && ax != az && ay == az){                                   
                              }else if(ax == ay && ax != az && ay != az){
                                  if(ax > az){return 0;}else if(ax < az){return 4;
                                  }else{System.out.println(" 0 a "+x+" - "+y+" - "+z+" Is not possible getQuadrant function.");}
                         //     }else if(ax != ay && ax == az && ay == az){                                   
                              }else if(ax != ay && ax == az && ay != az){ 
                                  if(ay > az){return 1;}else if(ay < az){return 4;
                                  }else{System.out.println(" 0 b "+x+" - "+y+" - "+z+" Is not possible getQuadrant function.");}
                              }else if(ax != ay && ax != az && ay == az){  
                                  if(ax > az){return 0;}else if(ax < az){return 4;
                                  }else{System.out.println(" 0 a "+x+" - "+y+" - "+z+" Is not possible getQuadrant function.");}
                        //      }else if(ax != ay && ax != az && ay != az){                                                                                                   
                              }else{System.out.println(" 0 z "+x+" - "+y+" - "+z+" Is not possible getQuadrant function.");}
                        }
               }else if(x >= 0 && y >= 0 && z < 0){  c1++;  // 1
                              if(ax > ay && ax > az && ay > az){return 0;    // 0
                        }else if(ax > ay && ax > az && ay < az){return 0;    // 1
                        }else if(ax > ay && ax < az && ay < az){return 5;    // 3
                        }else if(ax < ay && ax > az && ay > az){return 1;    // 4
                        }else if(ax < ay && ax < az && ay > az){return 1;    // 6
                        }else if(ax < ay && ax < az && ay < az){return 5;    // 7
                        }else{
                                    if(ax == ay && ax == az && ay == az){return 1;                                   
                              }else if(ax == ay && ax != az && ay != az){
                                  if(ax > az){return 1;}else if(ax < az){return 5;                                     
                                  }else{System.out.println(" 0 a "+x+" - "+y+" - "+z+" Is not possible getQuadrant function.");}                                  
                              }else if(ax != ay && ax == az && ay != az){ 
                                  if(ay > az){return 0;}else if(ay < az){return 5;                                     
                                  }else{System.out.println(" 0 b "+x+" - "+y+" - "+z+" Is not possible getQuadrant function.");}
                              }else if(ax != ay && ax != az && ay == az){  
                                  if(ax > az){return 1;}else if(ax < az){return 5;                                     
                                  }else{System.out.println(" 0 a "+x+" - "+y+" - "+z+" Is not possible getQuadrant function.");}                                                                                                
                              }else{System.out.println(" 0 z "+x+" - "+y+" - "+z+" Is not possible getQuadrant function.");}
                        }
               }else if(x >= 0 && y < 0 && z >= 0){  c2++;  // 2
                              if(ax > ay && ax > az && ay > az){return 0;    // 0
                        }else if(ax > ay && ax > az && ay < az){return 0;    // 1
                        }else if(ax > ay && ax < az && ay < az){return 4;    // 3
                        }else if(ax < ay && ax > az && ay > az){return 3;    // 4
                        }else if(ax < ay && ax < az && ay > az){return 3;    // 6
                        }else if(ax < ay && ax < az && ay < az){return 4;    // 7
                        }else{
                                    if(ax == ay && ax == az && ay == az){return 4;
                              }else if(ax == ay && ax != az && ay != az){
                                  if(ax > az){return 3;}else if(ax < az){return 4;
                                  }else{System.out.println(" 0 a "+x+" - "+y+" - "+z+" Is not possible getQuadrant function.");}
                              }else if(ax != ay && ax == az && ay != az){
                                  if(ay > az){return 0;}else if(ay < az){return 4;
                                  }else{System.out.println(" 0 b "+x+" - "+y+" - "+z+" Is not possible getQuadrant function.");}
                              }else if(ax != ay && ax != az && ay == az){
                                  if(ax > az){return 3;}else if(ax < az){return 4;
                                  }else{System.out.println(" 0 a "+x+" - "+y+" - "+z+" Is not possible getQuadrant function.");}
                              }else{System.out.println(" 0 z "+x+" - "+y+" - "+z+" Is not possible getQuadrant function.");}
                        }
               }else if(x >= 0 && y < 0 && z < 0){  c3++;  // 3
                              if(ax > ay && ax > az && ay > az){return 0;    // 0
                        }else if(ax > ay && ax > az && ay < az){return 0;    // 1
                        }else if(ax > ay && ax < az && ay < az){return 5;    // 3
                        }else if(ax < ay && ax > az && ay > az){return 3;    // 4
                        }else if(ax < ay && ax < az && ay > az){return 3;    // 6
                        }else if(ax < ay && ax < az && ay < az){return 5;    // 7
                        }else{                             
                                    if(ax == ay && ax == az && ay == az){return 5;                                   
                              }else if(ax == ay && ax != az && ay != az){
                                  if(ax > az){return 0;}else if(ax < az){return 5;                                     
                                  }else{System.out.println(" 0 a "+x+" - "+y+" - "+z+" Is not possible getQuadrant function.");}                                  
                              }else if(ax != ay && ax == az && ay != az){ 
                                  if(ay > az){return 3;}else if(ay < az){return 5;                                     
                                  }else{System.out.println(" 0 b "+x+" - "+y+" - "+z+" Is not possible getQuadrant function.");}
                              }else if(ax != ay && ax != az && ay == az){  
                                  if(ax > az){return 0;}else if(ax < az){return 5;                                     
                                  }else{System.out.println(" 0 a "+x+" - "+y+" - "+z+" Is not possible getQuadrant function.");}                                                                                                
                              }else{System.out.println(" 0 z "+x+" - "+y+" - "+z+" Is not possible getQuadrant function.");}
                        }
               }else if(x < 0 && y >= 0 && z >= 0){ c4++;   // 4
                              if(ax > ay && ax > az && ay > az){return 2;    // 0
                        }else if(ax > ay && ax > az && ay < az){return 2;    // 1
                        }else if(ax > ay && ax < az && ay < az){return 4;    // 3
                        }else if(ax < ay && ax > az && ay > az){return 1;    // 4
                        }else if(ax < ay && ax < az && ay > az){return 1;    // 6
                        }else if(ax < ay && ax < az && ay < az){return 4;    // 7
                        }else{
                                    if(ax == ay && ax == az && ay == az){return 4;
                              }else if(ax == ay && ax != az && ay != az){
                                  if(ax > az){return 1;}else if(ax < az){return 4;
                                  }else{System.out.println(" 0 a "+x+" - "+y+" - "+z+" Is not possible getQuadrant function.");}
                              }else if(ax != ay && ax == az && ay != az){
                                  if(ay > az){return 2;}else if(ay < az){return 4;
                                  }else{System.out.println(" 0 b "+x+" - "+y+" - "+z+" Is not possible getQuadrant function.");}
                              }else if(ax != ay && ax != az && ay == az){
                                  if(ax > az){return 1;}else if(ax < az){return 4;
                                  }else{System.out.println(" 0 a "+x+" - "+y+" - "+z+" Is not possible getQuadrant function.");}
                              }else{System.out.println(" 0 z "+x+" - "+y+" - "+z+" Is not possible getQuadrant function.");}
                        }
               }else if(x < 0 && y >= 0 && z < 0){  c5++;  // 5
                              if(ax > ay && ax > az && ay > az){return 2;    // 0
                        }else if(ax > ay && ax > az && ay < az){return 2;    // 1
                        }else if(ax > ay && ax < az && ay < az){return 5;    // 3
                        }else if(ax < ay && ax > az && ay > az){return 1;    // 4
                        }else if(ax < ay && ax < az && ay > az){return 1;    // 6
                        }else if(ax < ay && ax < az && ay < az){return 5;    // 7
                        }else{
                                    if(ax == ay && ax == az && ay == az){return 5;                                   
                              }else if(ax == ay && ax != az && ay != az){
                                  if(ax > az){return 2;}else if(ax < az){return 5;                                     
                                  }else{System.out.println(" 0 a "+x+" - "+y+" - "+z+" Is not possible getQuadrant function.");}                                  
                              }else if(ax != ay && ax == az && ay != az){ 
                                  if(ay > az){return 1;}else if(ay < az){return 5;                                     
                                  }else{System.out.println(" 0 b "+x+" - "+y+" - "+z+" Is not possible getQuadrant function.");}
                              }else if(ax != ay && ax != az && ay == az){  
                                  if(ax > az){return 2;}else if(ax < az){return 5;                                     
                                  }else{System.out.println(" 0 a "+x+" - "+y+" - "+z+" Is not possible getQuadrant function.");}                                                                                                
                              }else{System.out.println(" 0 z "+x+" - "+y+" - "+z+" Is not possible getQuadrant function.");}
                        }
               }else if(x < 0 && y < 0 && z >= 0){  c6++;  // 6
                              if(ax > ay && ax > az && ay > az){return 2;    // 0
                        }else if(ax > ay && ax > az && ay < az){return 2;    // 1
                        }else if(ax > ay && ax < az && ay < az){return 4;    // 3
                        }else if(ax < ay && ax > az && ay > az){return 3;    // 4
                        }else if(ax < ay && ax < az && ay > az){return 3;    // 6
                        }else if(ax < ay && ax < az && ay < az){return 4;    // 7
                        }else{
                                    if(ax == ay && ax == az && ay == az){return 3;                                   
                              }else if(ax == ay && ax != az && ay != az){
                                  if(ax > az){return 2;}else if(ax < az){return 4;                                     
                                  }else{System.out.println(" 0 a "+x+" - "+y+" - "+z+" Is not possible getQuadrant function.");}                                  
                              }else if(ax != ay && ax == az && ay != az){ 
                                  if(ay > az){return 3;}else if(ay < az){return 4;                                     
                                  }else{System.out.println(" 0 b "+x+" - "+y+" - "+z+" Is not possible getQuadrant function.");}
                              }else if(ax != ay && ax != az && ay == az){  
                                  if(ax > az){return 2;}else if(ax < az){return 4;                                     
                                  }else{System.out.println(" 0 a "+x+" - "+y+" - "+z+" Is not possible getQuadrant function.");}                                                                                                
                              }else{System.out.println(" 0 z "+x+" - "+y+" - "+z+" Is not possible getQuadrant function.");}
                        }
               }else if(x < 0 && y < 0 && z < 0){  c7++;  // 7
                              if(ax > ay && ax > az && ay > az){return 2;    // 0
                        }else if(ax > ay && ax > az && ay < az){return 2;    // 1
                        }else if(ax > ay && ax < az && ay < az){return 5;    // 3
                        }else if(ax < ay && ax > az && ay > az){return 3;    // 4
                        }else if(ax < ay && ax < az && ay > az){return 3;    // 6
                        }else if(ax < ay && ax < az && ay < az){return 5;    // 7
                        }else{
                                    if(ax == ay && ax == az && ay == az){return 2;                                   
                              }else if(ax == ay && ax != az && ay != az){
                                  if(ax > az){return 1;}else if(ax < az){return 4;                                     
                                  }else{System.out.println(" 0 a "+x+" - "+y+" - "+z+" Is not possible getQuadrant function.");}                                  
                              }else if(ax != ay && ax == az && ay != az){ 
                                  if(ay > az){return 2;}else if(ay < az){return 4;                                     
                                  }else{System.out.println(" 0 b "+x+" - "+y+" - "+z+" Is not possible getQuadrant function.");}
                              }else if(ax != ay && ax != az && ay == az){  
                                  if(ax > az){return 1;}else if(ax < az){return 4;                                     
                                  }else{System.out.println(" 0 a "+x+" - "+y+" - "+z+" Is not possible getQuadrant function.");}                                                                                                
                              }else{System.out.println(" 0 z "+x+" - "+y+" - "+z+" Is not possible getQuadrant function.");}
                        }
               }else{System.out.println(" e0 "+x+" - "+y+" - "+z+" Is not possible getQuadrant function.");}          
        return Q;
    }
    public int getQuadrant2(double x,double y, double z){
        int NoQuadrantAssaign = -1;
        if(y!=0){           
            if(x/y == -1){
              if(x>0){
                 if(x>z){return 3;}else{if(z>0){return 4;}else{return 5;}}  // q3 q4 q5 // o
              }else if(x<0){
                 if(x>z){return 1;}else{if(z>0){return 4;}else{return 5;}}  //q1 q4 q5  // g
              }
            }else if(x/y == 0){
              if(y>0){
                 if(y>z){return 1;}else{if(z>0){return 4;}else{return 5;}}  //q1 q4 q5  // e
              }else if(y<0){
                 if(y>z){return 3;}else{if(z>0){return 4;}else{return 5;}}  //q3 q4 q5  // m
              }
            }else if(x/y == 1){
              if(x>0){
                 if(x>z){return 0;}else{if(z>0){return 4;}else{return 5;}}  //q0 q4 q5  // c
              }else if(x<0){
                 if(x>z){return 2;}else{if(z>0){return 4;}else{return 5;}}  //q2 q4 q5  // k
              }
            }else{
               double absX = Math.abs(x);
               double absY = Math.abs(y);
               if(absX/absY>1){
                  if(x>0){
                     if(x>z){return 0;}else{if(z>0){return 4;}else{return 5;}} // q0 q4 q5 // b, p
                  }else{
                     if(x>z){return 2;}else{if(z>0){return 4;}else{return 5;}} // q2 14 15 // H, j
                  }
               }else{
                  if(y>0){
                     if(y>z){return 1;}else{if(z>0){return 4;}else{return 5;}} // q1 q4 q5 // f, d
                  }else{
                     if(y>z){return 3;}else{if(z>0){return 4;}else{return 5;}} // q3 q4 q5 // l, n
                  }
               }
            }
        }else{
           if(x==0){if(z > 0){return 4;}else if(z <= 0){return 5;}}
           if(z==0){if(x > 0){return 0;}else if(x <= 0){return 2;}} //  a, i        
        }
        return NoQuadrantAssaign;
    }
    public int getQuadrant3(double x, double y, double z){
        int Q = -1;
        double ax, ay, az;
        counter++;
        ax = Math.abs(Round((float)x,3));   ay = Math.abs(Round((float)y,3));   az = Math.abs(Round((float)z,3));
    //    System.out.println(counter+" A ("+ax+" "+ay+" "+az+")");
             if((ax != ay && ax != az && ay != az) && (ax != 0 && ay != 0 && az != 0)){
                     System.out.println(counter+" if ("+ax+" "+ay+" "+az+")");
                     if(x > 0 && y > 0 && z > 0){    // 0
                              if(ax > ay && ax > az && ay > az){return 0;    // 0  (3,2,1)
                        }else if(ax > ay && ax > az && ay < az){return 0;    // 1  (3,1,2)
                        }else if(ax > ay && ax < az && ay < az){return 4;    // 3  (,,,)
                        }else if(ax < ay && ax > az && ay > az){return 1;    // 4  (,,,)
                        }else if(ax < ay && ax < az && ay > az){return 1;    // 6  (,,,)
                        }else if(ax < ay && ax < az && ay < az){return 4;    // 7  (,,,)
                        }else{System.out.println(" 0 "+x+" - "+y+" - "+z+" Is not possible getQuadrant function.");}
               }else if(x > 0 && y > 0 && z < 0){    // 1
                              if(ax > ay && ax > az && ay > az){return 0;    // 0
                        }else if(ax > ay && ax > az && ay < az){return 0;    // 1
                        }else if(ax > ay && ax < az && ay < az){return 5;    // 3
                        }else if(ax < ay && ax > az && ay > az){return 1;    // 4
                        }else if(ax < ay && ax < az && ay > az){return 1;    // 6
                        }else if(ax < ay && ax < az && ay < az){return 5;    // 7
                        }else{System.out.println(" 1 "+x+" - "+y+" - "+z+" Is not possible getQuadrant function.");}
               }else if(x > 0 && y < 0 && z > 0){    // 2
                              if(ax > ay && ax > az && ay > az){return 0;    // 0
                        }else if(ax > ay && ax > az && ay < az){return 0;    // 1
                        }else if(ax > ay && ax < az && ay < az){return 4;    // 3
                        }else if(ax < ay && ax > az && ay > az){return 3;    // 4
                        }else if(ax < ay && ax < az && ay > az){return 3;    // 6
                        }else if(ax < ay && ax < az && ay < az){return 4;    // 7
                        }else{System.out.println(" 2 "+x+" - "+y+" - "+z+" Is not possible getQuadrant function.");}
               }else if(x > 0 && y < 0 && z < 0){    // 3
                              if(ax > ay && ax > az && ay > az){return 0;    // 0
                        }else if(ax > ay && ax > az && ay < az){return 0;    // 1
                        }else if(ax > ay && ax < az && ay < az){return 5;    // 3
                        }else if(ax < ay && ax > az && ay > az){return 3;    // 4
                        }else if(ax < ay && ax < az && ay > az){return 3;    // 6
                        }else if(ax < ay && ax < az && ay < az){return 5;    // 7
                        }else{System.out.println(" 3 "+x+" - "+y+" - "+z+" Is not possible getQuadrant function.");}
               }else if(x < 0 && y > 0 && z > 0){    // 4
                              if(ax > ay && ax > az && ay > az){return 2;    // 0
                        }else if(ax > ay && ax > az && ay < az){return 2;    // 1
                        }else if(ax > ay && ax < az && ay < az){return 4;    // 3
                        }else if(ax < ay && ax > az && ay > az){return 1;    // 4
                        }else if(ax < ay && ax < az && ay > az){return 1;    // 6
                        }else if(ax < ay && ax < az && ay < az){return 4;    // 7
                        }else{System.out.println(" 4 "+x+" - "+y+" - "+z+" Is not possible getQuadrant function.");}
               }else if(x < 0 && y > 0 && z < 0){    // 5
                              if(ax > ay && ax > az && ay > az){return 2;    // 0
                        }else if(ax > ay && ax > az && ay < az){return 2;    // 1
                        }else if(ax > ay && ax < az && ay < az){return 5;    // 3
                        }else if(ax < ay && ax > az && ay > az){return 1;    // 4
                        }else if(ax < ay && ax < az && ay > az){return 1;    // 6
                        }else if(ax < ay && ax < az && ay < az){return 5;    // 7
                        }else{System.out.println(" 5 "+x+" - "+y+" - "+z+" Is not possible getQuadrant function.");}
               }else if(x < 0 && y < 0 && z > 0){    // 6
                              if(ax > ay && ax > az && ay > az){return 2;    // 0
                        }else if(ax > ay && ax > az && ay < az){return 2;    // 1
                        }else if(ax > ay && ax < az && ay < az){return 4;    // 3
                        }else if(ax < ay && ax > az && ay > az){return 3;    // 4
                        }else if(ax < ay && ax < az && ay > az){return 3;    // 6
                        }else if(ax < ay && ax < az && ay < az){return 4;    // 7
                        }else{System.out.println(" 6 "+x+" - "+y+" - "+z+" IIs not possible getQuadrant function.");}
               }else if(x < 0 && y < 0 && z < 0){    // 7
                              if(ax > ay && ax > az && ay > az){return 2;    // 0
                        }else if(ax > ay && ax > az && ay < az){return 2;    // 1
                        }else if(ax > ay && ax < az && ay < az){return 5;    // 3
                        }else if(ax < ay && ax > az && ay > az){return 3;    // 4
                        }else if(ax < ay && ax < az && ay > az){return 3;    // 6
                        }else if(ax < ay && ax < az && ay < az){return 5;    // 7
                        }else{System.out.println(" 7 "+x+" - "+y+" - "+z+" Is not possible getQuadrant function.");}
               }else{System.out.println(" e0 "+x+" - "+y+" - "+z+" Is not possible getQuadrant function.");}
             }else{
                   System.out.println(counter+" M ("+ax+" "+ay+" "+az+")");

                           if(x == 0 && y == 0 && z != 0){
                              if(z > 0){return 4;}else{return 5;}
                     }else if(x == 0 && y != 0 && z == 0){
                              if(y > 0){return 1;}else{return 3;}
                     }else if(x == 0 && y != 0 && z != 0){
                              if(y > 0 && z > 0){
                                  if(ay > az){return 1;}else{return 4;}
                              }else if(y > 0 && z < 0){
                                  if(ay > az){return 1;}else{return 5;}
                              }else if(y < 0 && z > 0){
                                  if(ay > az){return 3;}else{return 4;}
                              }else if(y < 0 && z < 0){
                                  if(ay > az){return 3;}else{return 5;}
                              }
                     }else if(x != 0 && y == 0 && z == 0){
                              if(x > 0){return 0;}else{return 2;}
                     }else if(x != 0 && y == 0 && z != 0){
                              if(x > 0 && z > 0){
                                  if(ax > az){return 0;}else{return 4;}
                              }else if(x > 0 && z < 0){
                                  if(ax > az){return 0;}else{return 5;}
                              }else if(x < 0 && z > 0){
                                  if(ax > az){return 2;}else{return 4;}
                              }else if(x < 0 && z < 0){
                                  if(ax > az){return 2;}else{return 5;}
                              }
                     }else if(x != 0 && y != 0 && z == 0){
                              if(x > 0 && y > 0){
                                  if(ax > ay){return 0;}else{return 1;}
                              }else if(x > 0 && y < 0){
                                  if(ax > ay){return 0;}else{return 3;}
                              }else if(x < 0 && y > 0){
                                  if(ax > ay){return 2;}else{return 1;}
                              }else if(x < 0 && y < 0){
                                  if(ax > ay){return 2;}else{return 3;}
                              }
                     }else if(x != 0 && y != 0 && z != 0){
                         if(ax == ay && ax == az && ay == az){
                            if(ax > 0 && ay > 0 && az > 0){return 0;
                            }else if(ax > 0 && ay > 0 && az < 0){return 1;
                            }else if(ax > 0 && ay < 0 && az > 0){return 2;
                            }else if(ax > 0 && ay < 0 && az < 0){return 3;
                            }else if(ax < 0 && ay > 0 && az > 0){return 4;
                            }else if(ax < 0 && ay > 0 && az < 0){return 5;
                            }else if(ax < 0 && ay < 0 && az > 0){return 4;
                            }else if(ax < 0 && ay < 0 && az < 0){return 5;}
                         }else{System.out.println(" e2 "+ax+" | "+ay+" | "+az+" Is not possible getQuadrant function at this point.");}
                     }else{System.out.println(" e1 "+x+" - "+y+" - "+z+" Is not possible getQuadrant function.");}
            }
        return Q;
    }
    public void GenerateRandom(TheMatrix TM){
          double a1, a2, a3, f;
          double eAngx, eAngy, eAngz;
          String pattern = "##0.00000";
          DecimalFormat dimentionFormatter = new DecimalFormat(pattern);
          LOG.lout.println("n Ax Ay Az a1 a2 a3 q1 q2 q3 q4");
          for(int n = 0; n < TM.nMol; n++){
             eAngx = 2*Math.random()*Math.PI;
             eAngy = 2*Math.random()*Math.PI;
             eAngz = 2*Math.random()*Math.PI;
             LOG.lout.print(n+" "+dimentionFormatter.format(eAngx)+" "+dimentionFormatter.format(eAngy)+" "+dimentionFormatter.format(eAngz));
             a1 = 0.5 * eAngy;
             a2 = 0.5 * (eAngx - eAngz);
             a3 = 0.5 * (eAngx + eAngz);
             LOG.lout.print(" "+dimentionFormatter.format(a1)+" "+dimentionFormatter.format(a2)+" "+dimentionFormatter.format(a3));
             TM.q_u1[n] = Math.sin(a1) * Math.cos(a2);
             TM.q_u2[n] = Math.sin(a1) * Math.sin(a2);
             TM.q_u3[n] = Math.cos(a1) * Math.sin(a3);
             TM.q_u4[n] = Math.cos(a1) * Math.cos(a3);
             LOG.lout.println(" "+dimentionFormatter.format(TM.q_u1[n])+" "+dimentionFormatter.format(TM.q_u2[n])+
                             " "+dimentionFormatter.format(TM.q_u3[n])+" "+dimentionFormatter.format(TM.q_u4[n]));
          }
    }
    public void DefineMol(TheMatrix TM, String name){
        for(int j = 0; j < TM.sitesMol; j++){
            TM.rmx[j] = 0.0; TM.rmy[j] = 0.0; TM.rmz[j] = 0.0;
        }
        if(name.equalsIgnoreCase("tip3p")){
           TM.bCon = 120.3995;
           TM.ro = 3.15061;
           TM.ep = 0.1521;
           TM.rmz[0] = -0.0207;
           TM.rmz[1] = -0.0207;
           TM.rmy[2] = 0.2402;   TM.rmz[2] = 0.1653;
           TM.rmy[3] =-0.2402;   TM.rmz[3] = 0.1653;
           TM.typeF[0] = 1;   TM.typeF[1] = 2;    TM.typeF[2] = 3;    TM.typeF[3] = 3;
           TM.nAtoms = 3;     TM.nCharges = 3;
           // Only for tip4p it would have to be modified for other molecules or combination of molecules
           for(int j = 0; j < TM.nMol*TM.sitesMol; j=j+4){
               TM.atomType[j] = 1;
               TM.atomType[j+1] = 2;
               TM.atomType[j+2] = 3;
               TM.atomType[j+3] = 3;
           }
           for(int j = 0; j < TM.nMol; j++){TM.moleculeIndex[j] = (TM.sitesMol*3*j);}
        }else if(name.equalsIgnoreCase("tip4p")){
           TM.bCon = 183.6615;
           TM.ro = 3.15365;
           TM.ep = 0.1549;
           TM.rmz[0] =-0.02065;
           TM.rmz[1] = 0.02695;
           TM.rmy[2] = 0.2400;    TM.rmz[2] = 0.1652;
           TM.rmy[3] = -0.2400;   TM.rmz[3] = 0.1652;
           TM.typeF[0] = 1;   TM.typeF[1] = 2;    TM.typeF[2] = 3;    TM.typeF[3] = 3;
           TM.nAtoms = 3;      TM.nCharges = 3;
           // Only for tip4p it would have to be modified for other molecules or combination of molecules
           for(int j = 0; j < TM.nMol*TM.sitesMol; j=j+4){
               TM.atomType[j] = 1;
               TM.atomType[j+1] = 2;
               TM.atomType[j+2] = 3;
               TM.atomType[j+3] = 3;
           }
           for(int j = 0; j < TM.nMol; j++){TM.moleculeIndex[j] = (TM.sitesMol*3*j);}
        }else if(name.equalsIgnoreCase("tip5p")){
           TM.bCon = 39.4567;
           TM.ro = 3.1200;
           TM.ep = 0.1599;
           TM.rmz[0] = -0.0209;
           TM.rmy[1] = 0.2426;   TM.rmz[1] = 0.1669;
           TM.rmy[2] =-0.2426;   TM.rmz[2] = 0.1669;
           TM.rmx[3] = 0.1832;   TM.rmz[3] = -0.1504;
           TM.rmx[4] =-0.1832;   TM.rmz[4] = -0.1504;
           TM.typeF[0] = 1;   TM.typeF[1] = 3;    TM.typeF[2] = 3;
           TM.typeF[3] = 7;   TM.typeF[4] = 7;
           TM.nAtoms = 3;     TM.nCharges = 4;
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
           TM.bCon = 4.5224;
           TM.ro = 3.1000;
           TM.ep =0.0757;
           TM.rmx[0] = -0.0207;
           TM.rmy[1] = 0.263;   TM.rmz[1] = 0.1656;
           TM.rmy[2] =-0.263;   TM.rmz[2] = 0.1656;
           TM.rmx[3] = 0.2107;  TM.rmz[3] =-0.1697;
           TM.rmx[4] =-0.2107;  TM.rmz[4] =-0.1697;
           TM.typeF[0] = 1;   TM.typeF[1] = 3;    TM.typeF[2] = 3;    TM.typeF[3] = 7;   TM.typeF[4] = 7;
           TM.nAtoms = 3;     TM.nCharges = 4;
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
    public void GenSiteCoord(TheMatrix TM){
        double tx, ty, tz;
        for(int n = 0; n < TM.nMol; n++){
            for(int j = 0; j < TM.sitesMol; j++){
                tx = TM.rMatT[(n*9)+0]*TM.rmx[j] + TM.rMatT[(n*9)+3]*TM.rmy[j] + TM.rMatT[(n*9)+6]*TM.rmz[j];
                ty = TM.rMatT[(n*9)+1]*TM.rmx[j] + TM.rMatT[(n*9)+4]*TM.rmy[j] + TM.rMatT[(n*9)+7]*TM.rmz[j];
                tz = TM.rMatT[(n*9)+2]*TM.rmx[j] + TM.rMatT[(n*9)+5]*TM.rmy[j] + TM.rMatT[(n*9)+8]*TM.rmz[j];

                TM.rxs[n*TM.sitesMol + j] = TM.rx[n] + tx;
                TM.rys[n*TM.sitesMol + j] = TM.ry[n] + ty;
                TM.rzs[n*TM.sitesMol + j] = TM.rz[n] + tz;
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
}
