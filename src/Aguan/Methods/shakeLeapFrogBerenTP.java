package Aguan.Methods;

import Aguan.TheMatrix.TheMatrix;
import Aguan.files.restart;
import java.io.*;
import java.util.*;
import java.text.*;
/**
 * Time Step 0.00060241
 * @author  Noel Carrascal
 */
public class shakeLeapFrogBerenTP{
    private double[] initBox, finalBox;     int initBoxC, finalBoxC;
    private double[] initX, finalX;         int initXC, finalXC;
    private double[] initV, finalV;         int initVC, finalVC;
    private double[] initF, finalF;         int initFC, finalFC;
    private double[] Lxyz, Vxyz, Fxyz;
    private int nMol, sitesMol, stepLimit, NDIM, stepCount, nAtoms, nCharges, nPoints, vdwPoint;
    private int[] typeF, atomType, atomChrg, atomIndex,  moleculeIndex, moleculeAtoms;
    private double targetTemperature, rCut, deltaT, volume, velMag, timeNow, bCon, ro, ep, regionX, regionY, regionZ;
    private double eee, vdwe, vdwf;
    public shakeLeapFrogBerenTP(String[] args){
        nMol = 2;    sitesMol = 3;  // For tip3p
        // Some Parameters
        stepLimit = Integer.parseInt(args[5]);  NDIM = 3;    targetTemperature = 300;  //Kelvin
        rCut = 8; /**Angstoms*/   stepCount = 0;   deltaT = 2; /**femtoseconds*/ 
        velMag = Math.sqrt(NDIM*(1.0-1.0/nMol)*targetTemperature*0.5);
        timeNow = stepCount * deltaT;
        // Define parameters for tip3p
        typeF = new int[4];  //  atomType = new int[nMol*sitesMol];    
        bCon = 120.4953;     ro = 3.15061;          ep = 0.1520999;
        rCut = rCut*10/ro;   deltaT = deltaT*0.000625;
        typeF[0] = 1;   typeF[1] = 2;    typeF[2] = 3;    typeF[3] = 3;
        nAtoms = 3;     nCharges = 3;    nPoints = 3;     vdwPoint = 1;
        initBox = new double[(stepLimit+1)*NDIM*NDIM];            finalBox = new double[(stepLimit+1)*NDIM*NDIM];
          initX = new double[(stepLimit+1)*nMol*nAtoms*NDIM];       finalX = new double[(stepLimit+1)*nMol*nAtoms*NDIM];
          initV = new double[(stepLimit+1)*nMol*nAtoms*NDIM];       finalV = new double[(stepLimit+1)*nMol*nAtoms*NDIM];
          initF = new double[(stepLimit+1)*nMol*nAtoms*NDIM];       finalF = new double[(stepLimit+1)*nMol*nAtoms*NDIM];
        initBoxC = finalBoxC = initXC = finalXC = initVC = finalVC = initFC = finalFC = 0; 
        System.out.println(initBox.length+" "+initX.length+" "+initV.length+" "+initF.length);
        moleculeIndex = new int[nMol];  
        moleculeAtoms = new int[nMol];
        atomType = new int[nMol*nPoints]; /**1=O, 2,2=H */
        atomChrg = new int[nMol*nPoints];
        int molC = 0;
        for(int g = 0; g < atomType.length; g=g+nAtoms){
            atomType[g] = 1;     atomType[g+1] = atomType[g+2] = 2;
            atomChrg[g] = -2;    atomChrg[g+1] = atomChrg[g+2] = 1;
            moleculeIndex[molC] = g*NDIM; 
            moleculeAtoms[molC] = 3;  // For tip3p Only, might change for molecules with more atoms or vdw points.
            molC++;
        }
        Lxyz = new double[nMol*nPoints*NDIM];
        Vxyz = new double[nMol*nPoints*NDIM];
        Fxyz = new double[nMol*nPoints*NDIM];
        // Print All Arrays, and check visually
        int h;
        System.out.println("AtomType Array:");      for(h = 0; h < atomType.length; h++){ System.out.print(atomType[h]+", ");} System.out.println();
        System.out.println("AtomType Charge:");     for(h = 0; h < atomChrg.length; h++){ System.out.print(atomChrg[h]+", ");} System.out.println();
        System.out.println("moleculeIndex Array:"); for(h = 0; h < moleculeIndex.length; h++){ System.out.print(moleculeIndex[h]+", ");} System.out.println();
        System.out.println("moleculeAtoms Array:"); for(h = 0; h < moleculeAtoms.length; h++){ System.out.print(moleculeAtoms[h]+", ");} System.out.println();
        System.out.println("Lxyz Array:");          for(h = 0; h < Lxyz.length; h++){ System.out.print(Lxyz[h]+", ");} System.out.println();
        System.out.println("Vxyz Array:");          for(h = 0; h < Vxyz.length; h++){ System.out.print(Vxyz[h]+", ");} System.out.println();
        System.out.println("Fxyz Array:");          for(h = 0; h < Fxyz.length; h++){ System.out.print(Fxyz[h]+", ");} System.out.println();
        // TRR file generated with GROMACS is read
        readInput_trr(args[1]);
        // Assign Box in frame 0
        regionX = initBox[0]*(10/ro);     regionY = initBox[4]*(10/ro);      regionZ = initBox[8]*(10/ro); 
        // Assign frame 0 coords to temp arrays for processing
        for(int w = 0; w < (nMol*nAtoms*NDIM); w=w+3){
            Lxyz[w] = initX[w]*(10/ro);   Lxyz[w+1] = initX[w+1]*(10/ro);    Lxyz[w+2] = initX[w+2]*(10/ro);
            Vxyz[w] = initV[w];   Vxyz[w+1] = initV[w+1];    Vxyz[w+2] = initV[w+2];
        }
        System.out.println("Cxyz Box");
        System.out.println(regionX+" "+regionY+" "+regionZ);
        System.out.println("Lxyz Array:");  for(h = 0; h < Lxyz.length; h++){ System.out.print(Lxyz[h]+", ");} System.out.println();
        System.out.println("Vxyz Array:");  for(h = 0; h < Vxyz.length; h++){ System.out.print(Vxyz[h]+", ");} System.out.println();
        System.out.println("Fxyz Array:");  for(h = 0; h < Fxyz.length; h++){ System.out.print(Fxyz[h]+", ");} System.out.println();
        // Get forces at frame 0
        printFrame(0); 
        forces();
        printFrame(0);
        // Production Run 
        for(h = 1; h <= stepLimit; h++){

        }
     }
     public void forces(){
            int h,i,j,k,l;
            double dx, dy, dz, rr, rri, rri3, shiftx, shifty, shiftz;
            double eeeSum, vdweSum; 
            for(h = 0; h < (nMol*nAtoms*NDIM); h=h+3){
                Fxyz[h] = Fxyz[h+1] = Fxyz[h+2] = 0; 
            }
            eee = vdwe = eeeSum = 0; 
            for(i = 0; i < (nMol-1); i++){
                for(j = i+1; j < nMol; j++){
                    dx = Lxyz[moleculeIndex[i]] - Lxyz[moleculeIndex[j]];     
                    dy = Lxyz[moleculeIndex[i]+1] - Lxyz[moleculeIndex[j]+1];   
                    dz = Lxyz[moleculeIndex[i]+2] - Lxyz[moleculeIndex[j]+2];
                    rr = dx*dx + dy*dy + dz*dz;
                    shiftx = shifty = shiftz = 0;
                    if(dx >= 0.5*regionX)  shiftx = shiftx - regionX;
                    else if(dx < -0.5*regionX) shiftx = shiftx + regionX;
                    if(dy >= 0.5*regionY)  shifty = shifty - regionY;
                    else if(dy < -0.5*regionY) shifty = shifty + regionY;
                    if(dz >= 0.5*regionZ)  shiftz = shiftz - regionZ;
                    else if(dz < -0.5*regionZ) shiftz = shiftz + regionZ;
                    dx = dx + shiftx;  dy = dy + shifty;  dz = dz + shiftz;
                    System.out.printf("Shitxyz = %3.3f %3.3f %3.3f \n",shiftx, shifty, shiftz);
                    rr = dx*dx + dy*dy + dz*dz;  
                    if(rr < (rCut*rCut)){
                        for(k = moleculeIndex[i]; k < (moleculeIndex[i]+(moleculeAtoms[i]*NDIM)); k=k+NDIM){
                            for(l = moleculeIndex[j]; l < (moleculeIndex[j]+(moleculeAtoms[j]*NDIM)); l=l+NDIM){
                                dx = Lxyz[k] - Lxyz[l];     dy = Lxyz[k+1] - Lxyz[l+1];   dz = Lxyz[k+2] - Lxyz[l+2];
                                dx = dx + shiftx;           dy = dy + shifty;             dz = dz + shiftz; 
                                rr = dx*dx + dy*dy + dz*dz;
                                System.out.printf("k=%d %d,  l=%d %d -> %3.5f ",k,atomType[k/NDIM],l,atomType[l/NDIM],(Math.sqrt(rr)*ro));
                                rri = 1/rr;     
                                if((atomType[k/NDIM] == 1) && (atomType[l/NDIM] == 1)){
                                    rri3 = rri*rri*rri;
                                    vdwe = 4 * rri3 * (rri3 - 1);
                                    vdwf = 48 * rri3 * (rri3 - 0.5);
                                    eee = bCon*atomChrg[k/NDIM]*atomChrg[l/NDIM]*Math.sqrt(rri);
                                    eeeSum = eeeSum + eee;
                                    System.out.printf("  Energy vdwe = %3.10f    eee = %3.5f   eeeSum = %3.10f \n",(vdwe*ep*4.184),(eee*ep*4.184),(eeeSum*ep*4.184));
                                    finalF[k] = finalF[k] + (eee*dx*rri) + (vdwf*dx*rri);   finalF[k+1] = finalF[k+1] + (eee*dy*rri) + (vdwf*dy*rri);   finalF[k+2] = finalF[k+2] + (eee*dz*rri) + (vdwf*dz*rri);
                                    finalF[l] = finalF[l] - (eee*dx*rri) - (vdwf*dx*rri);   finalF[l+1] = finalF[l+1] - (eee*dy*rri) - (vdwf*dy*rri);   finalF[l+2] = finalF[l+2] - (eee*dz*rri) - (vdwf*dz*rri);
                                }else{
                                    eee = bCon*atomChrg[k/NDIM]*atomChrg[l/NDIM]*Math.sqrt(rri);
                                    eeeSum = eeeSum + eee;
                                    System.out.printf("  Energy eee = %3.5f eeeSum = %3.10f \n",(eee*ep*4.184),(eeeSum*ep*4.184));
                                    finalF[k] = finalF[k] + (eee*dx*rri);   finalF[k+1] = finalF[k+1] + (eee*dy*rri);   finalF[k+2] = finalF[k+2] + (eee*dz*rri);
                                    finalF[l] = finalF[l] - (eee*dx*rri);   finalF[l+1] = finalF[l+1] - (eee*dy*rri);   finalF[l+2] = finalF[l+2] - (eee*dz*rri);
                                }
                            }
                        } 
                    }
                }
            }
     }
     public void printFrame(int f){
        System.out.println("Frame "+f+"-------------------------------------");
        int w;
        System.out.println("Box(3x3)");
        for(w = (f*9); w < ((f*9)+9); w=w+3){
              System.out.printf("| %3.3f, %3.3f, %3.3f| - | %3.3f, %3.3f, %3.3f| =",initBox[w],initBox[w+1],initBox[w+2],finalBox[w],finalBox[w+1],finalBox[w+2]);
              System.out.printf("| %3.3f, %3.3f, %3.3f| \n",(initBox[w]-finalBox[w]),(initBox[w+1]-finalBox[w+1]),(initBox[w+2]-finalBox[w+2]));
        }
        System.out.println("xyz(6x3), vel(6x3), for(6x3)");
        int c = 1;
        for(w = (f*18); w < ((f*18)+18); w=w+3){
             System.out.printf("xyc[%d]| %3.3f, %3.3f, %3.3f| - | %3.3f, %3.3f, %3.3f| =",c,initX[w],initX[w+1],initX[w+2],finalX[w],finalX[w+1],finalX[w+2]);
             System.out.printf("| %3.3f, %3.3f, %3.3f| \n",(initX[w]-finalX[w]),(initX[w+1]-finalX[w+1]),(initX[w+2]-finalX[w+2]));
             c++;
        }
        c = 1;
        for(w = (f*18); w < ((f*18)+18); w=w+3){
             System.out.printf("vel[%d]| %3.3f, %3.3f, %3.3f| - | %3.3f, %3.3f, %3.3f| =",c,initV[w],initV[w+1],initV[w+2],finalV[w],finalV[w+1],finalV[w+2]);
             System.out.printf("| %3.3f, %3.3f, %3.3f| \n",(initV[w]-finalV[w]),(initV[w+1]-finalV[w+1]),(initV[w+2]-finalV[w+2]));
             c++;
        }
        c = 1;
        for(w = (f*18); w < ((f*18)+18); w=w+3){
             System.out.printf("FOR[%d]| %3.3f, %3.3f, %3.3f| - | %3.3f, %3.3f, %3.3f| =",c,(initF[w]),initF[w+1],initF[w+2],(finalF[w]*ep*4.184*ro),(finalF[w+1]*ep*4.184*ro),(finalF[w+2]*ep*4.184*ro));
             System.out.printf("| %3.3f, %3.3f, %3.3f| \n",(initF[w]-(finalF[w]*ep*4.184*ro)),(initF[w+1]-(finalF[w+1]*ep*4.184*ro)),(initF[w+2]-(finalF[w+2]*ep*4.184*ro)));
             c++;
        }
        System.out.println("Frame "+f+"-------------------------------------");
     }
    public void readInput_trr(String trr){
        try{
           FileReader FR =  new FileReader(trr);
           BufferedReader dataFile = new BufferedReader( FR );
           String line, temp, first;
           StringTokenizer token;
           int frameCount = 0;
           while( ( line = dataFile.readLine( ) ) != null ){
                    token = new StringTokenizer(line,":{}, ");
                    first = token.nextToken();
                    if(first.equals("natoms=")){
                      temp = token.nextToken();
                      temp = token.nextToken();  frameCount = Integer.parseInt(token.nextToken());
                    }else if(first.equals("box")){
                      for(int x = 0; x < 3; x++){
                          line = dataFile.readLine();
                          token = new StringTokenizer(line,":{}, ");
                          temp = token.nextToken();  temp = token.nextToken(); 
                          initBox[initBoxC] = Double.parseDouble(token.nextToken());  initBoxC++;
                          initBox[initBoxC] = Double.parseDouble(token.nextToken());  initBoxC++;
                          initBox[initBoxC] = Double.parseDouble(token.nextToken());  initBoxC++;
                          if(frameCount == 0){
                             finalBox[finalBoxC] = initBox[initBoxC-3];  finalBoxC++;
                             finalBox[finalBoxC] = initBox[initBoxC-2];  finalBoxC++;
                             finalBox[finalBoxC] = initBox[initBoxC-1];  finalBoxC++;
                          }
                      }
                    }else if(first.equals("x")){
                      for(int x = 0; x < (nAtoms*nMol); x++){
                          line = dataFile.readLine();
                          token = new StringTokenizer(line,":{}, ");
                          temp = token.nextToken();  temp = token.nextToken();
                          initX[initXC] = Double.parseDouble(token.nextToken());  initXC++;
                          initX[initXC] = Double.parseDouble(token.nextToken());  initXC++;
                          initX[initXC] = Double.parseDouble(token.nextToken());  initXC++;
                          if(frameCount == 0){
                             finalX[finalXC] = Lxyz[finalXC] = initX[initXC-3];  finalXC++;
                             finalX[finalXC] = Lxyz[finalXC] = initX[initXC-2];  finalXC++;
                             finalX[finalXC] = Lxyz[finalXC] = initX[initXC-1];  finalXC++;
                          }                   
                      }
                    }else if(first.equals("v")){
                      for(int x = 0; x < (nAtoms*nMol); x++){
                          line = dataFile.readLine();
                          token = new StringTokenizer(line,":{}, ");
                          temp = token.nextToken();  temp = token.nextToken();
                          initV[initVC] = Double.parseDouble(token.nextToken());  initVC++;
                          initV[initVC] = Double.parseDouble(token.nextToken());  initVC++;
                          initV[initVC] = Double.parseDouble(token.nextToken());  initVC++;
                          if(frameCount == 0){
                             finalV[finalVC] = Vxyz[finalVC] = initV[initVC-3];  finalVC++;
                             finalV[finalVC] = Vxyz[finalVC] = initV[initVC-2];  finalVC++;
                             finalV[finalVC] = Vxyz[finalVC] = initV[initVC-1];  finalVC++;
                          }      
                      }
                    }else if(first.equals("f")){
                      for(int x = 0; x < (nAtoms*nMol); x++){
                          line = dataFile.readLine();
                          token = new StringTokenizer(line,":{}, ");
                          temp = token.nextToken();  temp = token.nextToken();
                          initF[initFC] = Double.parseDouble(token.nextToken());  initFC++;
                          initF[initFC] = Double.parseDouble(token.nextToken());  initFC++;
                          initF[initFC] = Double.parseDouble(token.nextToken());  initFC++; 
                 //         if(frameCount == 0){
                 //            finalF[finalFC] = initF[initFC-3];  finalFC++;
                 //            finalF[finalFC] = initF[initFC-2];  finalFC++;
                 //            finalF[finalFC] = initF[initFC-1];  finalFC++;
                 //         }     
                      }
                   }
           }
        }catch( IOException e ){System.err.println( e );}
    }
/**
    public void LeapfrogStep(int part, TheMatrix  PrintManager PM){
          double[] mc, mt;
          double tx, ty, tz;
          mc = new double[9];          mt = new double[9];
          double temp = 0.5*deltaT;
          if(part == 1){
             for(int n = 0; n < nMol; n++){
                 rvx[n] = rvx[n] + temp * rax[n];
                 rvy[n] = rvy[n] + temp * ray[n];
                 rvz[n] = rvz[n] + temp * raz[n];
             }
             for(int n = 0; n < nMol; n++){
                 rx[n] = rx[n] + deltaT * rvx[n];
                 ry[n] = ry[n] + deltaT * rvy[n];
                 rz[n] = rz[n] + deltaT * rvz[n];
             }
          }else{
             for(int n = 0; n < nMol; n++){
                 rvx[n] = rvx[n] + temp * rax[n];
                 rvy[n] = rvy[n] + temp * ray[n];
                 rvz[n] = rvz[n] + temp * raz[n];
             }
          }
    }
    public void ComputeSiteForces(TheMatrix  PrintManager PM){
        double drx, dry, drz, shiftx, shifty, shiftz;
        double fcValEE, fcValRF, rr, sqrrr, rrCut, rrCut3, rri, rri3, uVal, uValVDW, uValEE, uValRF1, uValRF2;
        double fcValEEx, fcValEEy, fcValEEz, fcValRFx, fcValRFy, fcValRFz, virSumTemp;
        int ms1, ms2, typeSum;
        uValVDW = uValEE = uValRF1 = uValRF2 = fcValEE = fcValRF = 0;
        fcValEEx = fcValEEy = fcValEEz = fcValRFx = fcValRFy = fcValRFz = virSumTemp = 0;
        rrCut = rCut*rCut;
        rrCut3 = rrCut*rCut;
        for(int n = 0; n < nMol*sitesMol; n++){
              fxs[n] = 0.0; fys[n] = 0.0;  fzs[n] = 0.0;
        }
        uSum = uSumVDW = uSumEE = uSumRF1 = uSumRF2 = virSum = 0;
        for(int m1 = 0; m1 < nMol-1; m1++){
            for(int m2 = m1+1; m2 < nMol; m2++){
                drx = rx[m1] - rx[m2];
                dry = ry[m1] - ry[m2];
                drz = rz[m1] - rz[m2];
                shiftx = shifty = shiftz = 0.0;
                if(drx >= 0.5*regionX)  shiftx = shiftx - regionX;
                else if(drx < -0.5*regionX) shiftx = shiftx + regionX;
                if(dry >= 0.5*regionY)  shifty = shifty - regionY;
                else if(dry < -0.5*regionY) shifty = shifty + regionY;
                if(drz >= 0.5*regionZ)  shiftz = shiftz - regionZ;
                else if(drz < -0.5*regionZ) shiftz = shiftz + regionZ;
                drx = drx + shiftx;  dry = dry + shifty;  drz = drz + shiftz;
                rr = drx*drx + dry*dry + drz*drz;  
                sqrrr = Math.sqrt(rr);
                if(rr < rrCut){
                    ms1 = m1 * sitesMol;
                    ms2 = m2 * sitesMol;
                    //System.out.printf("%d m1(%d)-m2(%d) ",stepCount,m1,m2);
                    for(int j1 = 0; j1 < sitesMol ; j1++){
                        for(int j2 = 0; j2 < sitesMol ; j2++){
                            uValVDW = uValEE = uValRF1 = uValRF2 = fcValEE = fcValRF = 0;
                            typeSum =  typeF[j1] + typeF[j2];
                            if(typeF[j1] == typeF[j2] || typeSum == 5 || typeSum == 10){
                                drx = rxs[ms1+j1] - rxs[ms2+j2];
                                dry = rys[ms1+j1] - rys[ms2+j2];
                                drz = rzs[ms1+j1] - rzs[ms2+j2];
                                drx = drx + shiftx;  dry = dry + shifty;  drz = drz + shiftz;
                                rr = drx*drx + dry*dry + drz*drz;                                
                                rri = 1.0/rr;                                
                                switch(typeSum){
                                    case 2:
                                        rri3 = rri*rri*rri;
                                        uValVDW = 4 * rri3 * (rri3 - 1);
                                        virSumTemp = 48 * rri3 * (rri3 - 0.5);
                                        fcValEE = virSumTemp * rri;
                                        break;
                                    case 4:
                                        uValEE = 4 * bCon * Math.sqrt(rri);
                                        uValRF1 = 4 * bCon * rr/rrCut3;
                                        uValRF2 = bCon * (2*nMol*((-1*rxs[ms1+j1]*rxs[ms2+j2])+(-1*rys[ms1+j1]*rys[ms2+j2])+(-1*rzs[ms1+j1]*rzs[ms2+j2])))/rrCut3;                                        
                                        virSumTemp = uValEE + uValRF1;
                                        fcValEE = (uValEE * rri);
                                        fcValRF = (uValRF1/rrCut3);
                                        break;
                                    case 5:
                                        uValEE = -2 * bCon * Math.sqrt(rri);
                                        uValRF1 = -2 * bCon * rr/rrCut3;
                                        uValRF2 = bCon * (-nMol*((-1*rxs[ms1+j1]*rxs[ms2+j2])+(-1*rys[ms1+j1]*rys[ms2+j2])+(-1*rzs[ms1+j1]*rzs[ms2+j2])))/rrCut3;                                        
                                        virSumTemp = uValEE + uValRF1;
                                        fcValEE = (uValEE * rri);
                                        fcValRF = (uValRF1/rrCut3);
                                        break;
                                    case 6:                                       
                                        uValEE = bCon * Math.sqrt(rri);
                                        uValRF1 = bCon * rr/rrCut3;
                                        uValRF2 = bCon * (0.5*nMol*((-1*rxs[ms1+j1]*rxs[ms2+j2])+(-1*rys[ms1+j1]*rys[ms2+j2])+(-1*rzs[ms1+j1]*rzs[ms2+j2])))/rrCut3;
                                        virSumTemp = uValEE + uValRF1;
                                        fcValEE = (uValEE * rri);
                                        fcValRF = (uValRF1/rrCut3);
                                        break;
                                    case 10:                                        
                                        uValEE = -1 * bCon * Math.sqrt(rri);
                                        uValRF1 = -1 * bCon * rr/rrCut3;
                                        uValRF2 = bCon * (-0.5*nMol*((-1*rxs[ms1+j1]*rxs[ms2+j2])+(-1*rys[ms1+j1]*rys[ms2+j2])+(-1*rzs[ms1+j1]*rzs[ms2+j2])))/rrCut3;
                                        virSumTemp = uValEE + uValRF1;
                                        fcValEE = (uValEE * rri);
                                        fcValRF = (uValRF1/rrCut3);
                                        break;
                                    case 14:                                        
                                        uValEE = bCon * Math.sqrt(rri);
                                        uValRF1 = bCon * rr/rrCut3;
                                        uValRF2 = bCon * (0.5*nMol*((-1*rxs[ms1+j1]*rxs[ms2+j2])+(-1*rys[ms1+j1]*rys[ms2+j2])+(-1*rzs[ms1+j1]*rzs[ms2+j2])))/rrCut3;
                                        virSumTemp = uValEE + uValRF1;
                                        fcValEE = (uValEE * rri);
                                        fcValRF = (uValRF1/rrCut3);
                                        break;
                                }
                                fcValEEx = fcValEE*drx;  fcValEEy = fcValEE*dry;     fcValEEz = fcValEE*drz;
                                fcValRFx = fcValRF*drx;  fcValRFy = fcValRF*dry;     fcValRFz = fcValRF*drz;
                                virSum = virSum - (virSumTemp);
                                fxs[ms1+j1] = fxs[ms1+j1] + fcValEEx + fcValRFx;
                                fys[ms1+j1] = fys[ms1+j1] + fcValEEy + fcValRFy;
                                if(typeF[j1] == 1){
                                   fzs[ms1+j1] = fzs[ms1+j1] + fcValEEz + fcValRFz;
                                }else if(typeF[j1] == 2){
                                   fzs[ms1+j1] = fzs[ms1+j1] + fcValEEz + fcValRFz + (electricField*(-2));
                                }else if(typeF[j1] == 3){
                                   fzs[ms1+j1] = fzs[ms1+j1] + fcValEEz + fcValRFz + (electricField);
                                }else if(typeF[j1] == 7){
                                   fzs[ms1+j1] = fzs[ms1+j1] + fcValEEz + fcValRFz + (electricField*(-1));
                                }
                                
                                fxs[ms2+j2] = fxs[ms2+j2] - fcValEEx - fcValRFx;
                                fys[ms2+j2] = fys[ms2+j2] - fcValEEy - fcValRFy;
                                if(typeF[j2] == 1){
                                   fzs[ms2+j2] = fzs[ms2+j2] - fcValEEz - fcValRFz;
                                }else if(typeF[j2] == 2){
                                   fzs[ms2+j2] = fzs[ms2+j2] - fcValEEz  - fcValRFz + (electricField*(-2));
                                }else if(typeF[j2] == 3){
                                   fzs[ms2+j2] = fzs[ms2+j2] - fcValEEz  - fcValRFz + (electricField);
                                }else if(typeF[j2] == 7){
                                   fzs[ms2+j2] = fzs[ms2+j2] - fcValEEz  - fcValRFz + (electricField*(-1));
                                }
                                
                                uSumVDW += uValVDW;
                                uSumEE += uValEE;
                                uSumRF1 += uValRF1;
                                uSumRF2 += uValRF2;
                                uSum += (uValVDW + uValEE + uValRF1 + uValRF2);                                     
                            }
                        }
                    }                  
                }
            }
        }//A one molecule test could be commented out for multiple molecules, longer simulations
         //to improve performance.
        if(nMol == 1){
           for(int j1 = 0; j1 < sitesMol ; j1++){
               
           }
        }
        //PM.lout.printf("Sum ENER:\n");
        //PM.lout.printf("%3.3f %3.3f %3.3f %3.3f \n",uSumVDW,uSumEE,uSumRF1,uSumRF2);
    }
    public void ApplyBoundaryCond(TheMatrix {
         for(int n = 0; n < nMol; n++){
                if( rx[n] >= 0.5 * regionX){rx[n] = rx[n] - regionX;}
           else if( rx[n] < -0.5 * regionX){rx[n] = rx[n] + regionX;}
                if( ry[n] >= 0.5 * regionY){ry[n] = ry[n] - regionY;}
           else if( ry[n] < -0.5 * regionY){ry[n] = ry[n] + regionY;}
                if( rz[n] >= 0.5 * regionZ){rz[n] = rz[n] - regionZ;}
           else if( rz[n] < -0.5 * regionZ){rz[n] = rz[n] + regionZ;}
         }
    }
    public void AdjustTemp(TheMatrix  PrintManager PM){
        double wvBx, wvBy, wvBz;
        double vFac, vFacT, vFacR;
        //   vFacT = Math.sqrt(1+(deltaT/tau)*((targetTemperature/(instantTemperature*tranFrac)) - 1));
        //   vFacR = Math.sqrt(1+(deltaT/tau)*((targetTemperature/(instantTemperature*rotFrac)) - 1));
        if(ThermostatType == 1){  // According with the Equipartition, No Berensen.
               vFacT = velMag/Math.sqrt(vvSum/(nMol/2));
               vvSum = 0.0;
               for(int n = 0; n < nMol; n++){
                   rvx[n] *= vFacT;     rvy[n] *= vFacT;    rvz[n] *= vFacT;
                   vvSum += rvx[n]*rvx[n] + rvy[n]*rvy[n] + rvz[n]*rvz[n];
               }
               translationalTemperature = 2*(vvSum/nMol)/3.0;
               trzKinEnergyVal = 0.5 * vvSum/nMol;
               vFacR = velMag/Math.sqrt(vvrSum/(nMol/2));
               vvrSum = 0.0;
               for(int n = 0; n < nMol; n++){
                   wvx[n] *= vFacR;     wvy[n] *= vFacR;    wvz[n] *= vFacR;
                   wvBx = rMatT[(n*9)+0]*wvx[n] + rMatT[(n*9)+1]*wvy[n] + rMatT[(n*9)+2]*wvz[n];
                   wvBy = rMatT[(n*9)+3]*wvx[n] + rMatT[(n*9)+4]*wvy[n] + rMatT[(n*9)+5]*wvz[n];
                   wvBz = rMatT[(n*9)+6]*wvx[n] + rMatT[(n*9)+7]*wvy[n] + rMatT[(n*9)+8]*wvz[n];
                   vvrSum += mInertX*wvBx*wvBx + mInertY*wvBy*wvBy + mInertZ*wvBz*wvBz;
               }
               rotationalTemperature = 2*(vvrSum/nMol)/3.0;
               rotKinEnergyVal = 0.5 * vvrSum/nMol;
        }else if(ThermostatType == 2){
               vFacT = Math.sqrt(1+(tau*(((0.25*targetTemperature)/translationalTemperature)-1)));
               vvSum = 0.0;
               for(int n = 0; n < nMol; n++){
                   rvx[n] *= vFacT;     rvy[n] *= vFacT;    rvz[n] *= vFacT;
                   vvSum += rvx[n]*rvx[n] + rvy[n]*rvy[n] + rvz[n]*rvz[n];
               }
               translationalTemperature = 2*(vvSum/nMol)/3.0;
               trzKinEnergyVal = 0.5 * vvSum/nMol;
               vFacR = Math.sqrt(1+(tau*(((0.25*targetTemperature)/rotationalTemperature)-1)));
               vvrSum = 0.0;
               for(int n = 0; n < nMol; n++){
                   wvx[n] *= vFacR;     wvy[n] *= vFacR;    wvz[n] *= vFacR;
                   wvBx = rMatT[(n*9)+0]*wvx[n] + rMatT[(n*9)+1]*wvy[n] + rMatT[(n*9)+2]*wvz[n];
                   wvBy = rMatT[(n*9)+3]*wvx[n] + rMatT[(n*9)+4]*wvy[n] + rMatT[(n*9)+5]*wvz[n];
                   wvBz = rMatT[(n*9)+6]*wvx[n] + rMatT[(n*9)+7]*wvy[n] + rMatT[(n*9)+8]*wvz[n];
                   vvrSum += mInertX*wvBx*wvBx + mInertY*wvBy*wvBy + mInertZ*wvBz*wvBz;
               }
               rotationalTemperature = 2*(vvrSum/nMol)/3.0;
               rotKinEnergyVal = 0.5 * vvrSum/nMol;
        }else if(ThermostatType == 3){
               vFac = Math.sqrt(1+(tau*(((0.5*targetTemperature)/instantTemperature)-1)));
               vvSum = 0.0;
               for(int n = 0; n < nMol; n++){
                   rvx[n] *= vFac;     rvy[n] *= vFac;    rvz[n] *= vFac;
                   vvSum += rvx[n]*rvx[n] + rvy[n]*rvy[n] + rvz[n]*rvz[n];
               }
               translationalTemperature = 2*(vvSum/nMol)/3.0;
               trzKinEnergyVal = 0.5 * vvSum/nMol;
               vvrSum = 0.0;
               for(int n = 0; n < nMol; n++){
                   wvx[n] *= vFac;     wvy[n] *= vFac;    wvz[n] *= vFac;
                   wvBx = rMatT[(n*9)+0]*wvx[n] + rMatT[(n*9)+1]*wvy[n] + rMatT[(n*9)+2]*wvz[n];
                   wvBy = rMatT[(n*9)+3]*wvx[n] + rMatT[(n*9)+4]*wvy[n] + rMatT[(n*9)+5]*wvz[n];
                   wvBz = rMatT[(n*9)+6]*wvx[n] + rMatT[(n*9)+7]*wvy[n] + rMatT[(n*9)+8]*wvz[n];
                   vvrSum += mInertX*wvBx*wvBx + mInertY*wvBy*wvBy + mInertZ*wvBz*wvBz;
               }
               rotationalTemperature = 2*(vvrSum/nMol)/3.0;
               rotKinEnergyVal = 0.5 * vvrSum/nMol;
        }else{
              vvSum = 0.0;
              for(int n = 0; n < nMol; n++){
                  rvx[n] *= rvxf[n];     rvy[n] *= rvyf[n];    rvz[n] *= rvzf[n];
                  vvSum += rvx[n]*rvx[n] + rvy[n]*rvy[n] + rvz[n]*rvz[n];
              }
              translationalTemperature = (vvSum/nMol)/3.0;
              trzKinEnergyVal = 0.5 * vvSum/nMol;

              vvrSum = 0.0;
              for(int n = 0; n < nMol; n++){
                  wvx[n] *= rrvxf[n];     wvy[n] *= rrvyf[n];    wvz[n] *= rrvzf[n];
                  wvBx = rMatT[(n*9)+0]*wvx[n] + rMatT[(n*9)+1]*wvy[n] + rMatT[(n*9)+2]*wvz[n];
                  wvBy = rMatT[(n*9)+3]*wvx[n] + rMatT[(n*9)+4]*wvy[n] + rMatT[(n*9)+5]*wvz[n];
                  wvBz = rMatT[(n*9)+6]*wvx[n] + rMatT[(n*9)+7]*wvy[n] + rMatT[(n*9)+8]*wvz[n];
                  vvrSum += mInertX*wvBx*wvBx + mInertY*wvBy*wvBy + mInertZ*wvBz*wvBz;
              }
              rotationalTemperature = (vvrSum/nMol)/3.0;
              rotKinEnergyVal = 0.5 * vvrSum/nMol;
        }
        kinEnergyVal = rotKinEnergyVal + trzKinEnergyVal;
        instantTemperature = translationalTemperature + rotationalTemperature;
        totEnergyVal = kinEnergyVal + potEnergyVal;
    }
    public void EvalProps(TheMatrix  PrintManager PM){
     //    PM.lout.printf("Eval Prop (IN) \n");
         double wvBx, wvBy, wvBz;
         wvBx = wvBy = wvBz = 0;
         vSumX = vSumY = vSumZ = 0.0;
         vvSum = vvrSum = 0.0;
         double TKE = 0.0; double WKE = 0.0;
         for(int n = 0; n < nMol; n++){
         //     vSumX = vSumX + rvx[n];   // This comented lines could be used for
         //     vSumY = vSumY + rvy[n];   // correcting drifting in MD and possibly
         //     vSumZ = vSumZ + rvz[n];   // to avoid the 'flying ice problem'
             vvSum += ((rvx[n]*rvx[n]) + (rvy[n]*rvy[n]) + (rvz[n]*rvz[n]));
         }
         translationalTemperature = (vvSum/nMol)/3.0;
         instantTemperature = translationalTemperature + rotationalTemperature;
         trzKinEnergyVal = 0.5 * vvSum/nMol;
         kinEnergyVal = trzKinEnergyVal;         
         //pressure = density*(vvSum + vvrSum + virSum)/(nMol*NDIM);
         pressure = ((vvSum + vvrSum) - virSum)/(NDIM*nMol*volume);
         potEnergyVdwVal = uSumVDW/nMol;
         potEnergyEeVal = uSumEE/nMol;
         potEnergyRf1Val = uSumRF1/nMol;
         potEnergyRf2Val = uSumRF2/nMol;         
         potEnergyVal = potEnergyVdwVal + potEnergyEeVal + potEnergyRf1Val + potEnergyRf2Val + potEnergyEfVal;
         totEnergyVal = kinEnergyVal + potEnergyVal;
         tranFrac = 0.5;
         if(ThermostatType == 4){
             tranFrac = translationalTemperature/instantTemperature;
             for(int n = 0; n < nMol; n++){
                 rvxf[n] = rvx[n]/vvSum;
                 rvyf[n] = rvy[n]/vvSum;
                 rvzf[n] = rvz[n]/vvSum;
             }
         }
    }
    public void AccumProps(TheMatrix  int icode){
        if(icode == 0){
            totEnergySum = kinEnergySum = potEnergySum = 0.0;
            potEnergyVdwSum = potEnergyEeSum = potEnergyRf1Sum = potEnergyRf2Sum = potEnergyEfSum = 0;
            trzKinEnergySum = rotKinEnergySum = 0.0;
            aveTemperature = pressure = pressureSum = 0;
            rotTempAve = trzTempAve = 0;
        }else if(icode == 1){
            totEnergySum = totEnergySum + totEnergyVal;
            potEnergySum = potEnergySum + potEnergyVal;
            potEnergyVdwSum = potEnergyVdwSum + potEnergyVdwVal;
            potEnergyEeSum = potEnergyEeSum + potEnergyEeVal;
            potEnergyRf1Sum = potEnergyRf1Sum + potEnergyRf1Val;
            potEnergyRf2Sum = potEnergyRf2Sum + potEnergyRf2Val;
            potEnergyEfSum = potEnergyEfSum + potEnergyEfVal;
            kinEnergySum = kinEnergySum + kinEnergyVal;
            trzKinEnergySum = trzKinEnergySum + trzKinEnergyVal;
            rotKinEnergySum = rotKinEnergySum + rotKinEnergyVal;
            aveTemperature = aveTemperature + instantTemperature;
            rotTempAve = rotTempAve + rotationalTemperature;
            trzTempAve = trzTempAve + translationalTemperature;
            pressureSum = pressureSum + pressure;
        }else if(icode == 2){
            totEnergySum /= stepAvg;
            potEnergySum /= stepAvg;
            potEnergyVdwSum /= stepAvg;
            potEnergyEeSum /= stepAvg;
            potEnergyRf1Sum /= stepAvg;
            potEnergyRf2Sum /= stepAvg;
            potEnergyEfSum /= stepAvg;
            kinEnergySum /= stepAvg;
            trzKinEnergySum /= stepAvg;
            rotKinEnergySum /= stepAvg;
            aveTemperature /= stepAvg;
            rotTempAve /= stepAvg;
            trzTempAve /= stepAvg;
            pressureSum /= stepAvg;
        }
    }
*/
}
