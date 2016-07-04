/*
* To change this template, choose Tools | Templates
* and open the template in the editor.
*/

package Aguan.files;

import Aguan.TheMatrix.TheMatrix;
import java.io.*;
import java.util.*;
import Aguan.files.*;
import Aguan.parameters;
/**
 * @author Noel
 */
public class restart extends file {
    public PrintWriter pout;
    private TheMatrix TM;
    private double pwer;
    private Random generator;
    private String file_out; 
    private psf PSF;
    private crd CRD;
    private gro GRO;
/**    public restart(String in, String out){
         super(in);
         super.file_type = "restart";
         file_out = out;
         super.openOutputWriter(file_out);
         pout = super.PW;
    }*/
    public restart(String in){
         super(in);
         super.file_type = "restart";
    }
    public restart(String[] args){
         super();
         super.file_type = "restart";
         TM = new TheMatrix();
      if(args.length == 1){
         System.out.println("Follow -c by the following restart commands:");
         command_description();
      }else{
         if(args[1].equals("-random")){
               super.openOutputWriter(args[2]+"_0.rst");
               pout = super.PW;
               TM.restartCount = 0;
               TM.nMol = Integer.parseInt(args[3]);         // Number of water molecules.
               parameters.DefinePar(TM,args[4]);            // Water model.
               TM.initMatrix();
               parameters.DefineMol(TM,args[4]);     
               TM.regionX = Double.parseDouble(args[5])/TM.ro;    // Box X
               TM.regionY = Double.parseDouble(args[6])/TM.ro;    // Box Y
               TM.regionZ = Double.parseDouble(args[7])/TM.ro;    // Box Z
               RandomGenerate(TM, Long.parseLong(args[8]), args[1]); // args[8] is the seed to reproduce random number generation.
               writeRestartIn(TM);
        }else if(args[1].equals("-crystal")){
               super.openOutputWriter(args[2]+"_0.rst");
               pout = super.PW;
               TM.restartCount = 0;
               TM.nMol = Integer.parseInt(args[3]);
               parameters.DefinePar(TM,args[4]);
               TM.initMatrix();
               parameters.DefineMol(TM,args[4]);
               TM.density = Double.parseDouble(args[5])*(TM.density_convFact);
               TM.initUcellX = Integer.parseInt(args[6]);
               TM.initUcellY = Integer.parseInt(args[7]);
               TM.initUcellZ = Integer.parseInt(args[8]);
               TM.regionX = TM.initUcellX/Math.cbrt(TM.density);
               TM.regionY = TM.initUcellY/Math.cbrt(TM.density);
               TM.regionZ = TM.initUcellZ/Math.cbrt(TM.density);
               RandomGenerate(TM, Long.parseLong(args[9]), args[1]);  // args[9] is the seed to reproduce random number generation.
               writeRestartIn(TM);
         }else if(args[1].equals("-psf")){
               TM.restartIn = args[2];
               ReadRestartHeader(TM);
               TM.initMatrix();
               parameters.DefineMol(TM,TM.molType);
               ReadRestartFile(TM);
               PSF = new psf(args[2].substring(0,args[2].length()-4)+".psf");
               PSF.createPSFFile(TM);
         }else if(args[1].equals("-crd")){
               TM.restartIn = args[2];
               ReadRestartHeader(TM);
               TM.initMatrix();
               parameters.DefineMol(TM,TM.molType);
               ReadRestartFile(TM);
               GenSiteCoord(TM);
               CRD = new crd(args[2].substring(0,args[2].length()-4)+".crd");
               CRD.createCRDFile(TM);
         }else if(args[1].equals("-gro")){
               TM.restartIn = args[2];
               ReadRestartHeader(TM);
               TM.initMatrix();
               parameters.DefineMol(TM,TM.molType);
               ReadRestartFile(TM);
               GenSiteCoord(TM);
               GRO = new gro(args[2].substring(0,args[2].length()-4)+".gro");
               GRO.createGROFile(TM);
         }else if(args[1].equals("-rescale")){
               //ReadRestartFile(TM);
               //TM.oldDensity = TM.density;
               //TM.set_molType_sitesMol(args[4]);
               //TM.density = Double.parseDouble(args[5]);
               //Rescale(TM);
         }
      }
    }
    public void setOutputFile(String out){
         file_out = out;
         super.openOutputWriter(file_out);
         pout = super.PW;
    }
    public void Rescale(TheMatrix TM){
         double volume2, regionX2, regionY2, regionZ2;
         TM.regionX = (1/Math.pow(TM.density,pwer))*TM.initUcellX;
         TM.regionY = (1/Math.pow(TM.density,pwer))*TM.initUcellY;
         TM.regionZ = (1/Math.pow(TM.density,pwer))*TM.initUcellZ;

         regionX2 = (1/Math.pow(TM.oldDensity,pwer))*TM.initUcellX;
         regionY2 = (1/Math.pow(TM.oldDensity,pwer))*TM.initUcellY;
         regionZ2 = (1/Math.pow(TM.oldDensity,pwer))*TM.initUcellZ;

         TM.volume = TM.regionX*TM.regionY*TM.regionZ;
         volume2 = regionX2*regionY2*regionZ2;
         System.out.printf("Adjusted Volume = %3.5f density = %3.4f \n",TM.volume,Math.pow(TM.density,pwer));
         System.out.printf("Starting Volume = %3.5f density = %3.4f \n",volume2,Math.pow(TM.oldDensity,pwer));
         System.out.println("Box size was changed. ");
         System.out.printf("from RegionesXYZ = %3.3f %3.3f %3.3f \n",regionX2,regionY2,regionZ2);
         System.out.printf("to RegionesXYZ2 = %3.3f %3.3f %3.3f \n",TM.regionX,TM.regionY,TM.regionZ);
         // This means thet the molecules are rescaled down from the old density
         // and then scaled back up to the new density.
         // There should be a flag for this. One may want to do either experiment.
         // Also, one might want to rescale density by moving only one dimension,
         // as in the case of a piston.
         for(int n=0; n < TM.nMol; n++){
             TM.rx[n] /= regionX2;
             TM.ry[n] /= regionY2;
             TM.rz[n] /= regionZ2;
         }
         for(int n=0; n < TM.nMol; n++){
             TM.rx[n] *= TM.regionX;
             TM.ry[n] *= TM.regionY;
             TM.rz[n] *= TM.regionZ;
         }
    }
    public void writeRestartIn(TheMatrix TM){
            for(int i = 0; i < 12; i++){
                if(i == 0){
                    pout.println("STEP "+TM.stepCount+" "+TM.restartCount);
                }else if(i == 1){
                    pout.println("MOLECULES "+TM.nMol+" "+TM.molType);
                }else if(i == 2){
                    pout.println("BOX "+(TM.regionX*TM.ro)+" "+(TM.regionY*TM.ro)+" "+(TM.regionZ*TM.ro));
                }else if(i == 3){
                    pout.println("COORDS n= "+TM.nMol+" x 3 X(0),Y(0),Z(0)...X(n-1),Y(n-1),Z(n-1)");
                    for(int j = 0; j < TM.nMol; j++){
                        pout.println(j+" "+(TM.rx[j]*TM.ro)+" "+(TM.ry[j]*TM.ro)+" "+(TM.rz[j]*TM.ro));
                    }
                }else if(i == 4){
                    pout.println("VELS n= "+TM.nMol+" x 3 X(0),Y(0),Z(0)...X(n-1),Y(n-1),Z(n-1)");
                    for(int j = 0; j < TM.nMol; j++){
                       pout.println(j+" "+TM.rvx[j]+" "+TM.rvy[j]+" "+TM.rvz[j]);
                    }
                }else if(i == 5){
                    pout.println("ACCELS n= "+TM.nMol+" x 3 X(0),Y(0),Z(0)...X(n-1),Y(n-1),Z(n-1)");
                    for(int j = 0; j < TM.nMol; j++){
                       pout.println(j+" "+TM.rax[j]+" "+TM.ray[j]+" "+TM.raz[j]);
                    }
                }else if(i == 6){
                    pout.println("ANGCOORDS n= "+TM.nMol+" x 4 q1(0),q2(0),q3(0),q(4)....q1(n-1),q2(n-1),q3(n-1),q4(n-1)");
                    for(int j = 0; j < TM.nMol; j++){
                       pout.print(j+" "+TM.rMatT[j*9]+" "+TM.rMatT[j*9+1]+" "+TM.rMatT[j*9+2]);
                       pout.print(" "+TM.rMatT[j*9+3]+" "+TM.rMatT[j*9+4]+" "+TM.rMatT[j*9+5]);
                       pout.println(" "+TM.rMatT[j*9+6]+" "+TM.rMatT[j*9+7]+" "+TM.rMatT[j*9+8]);
                    }
                }else if(i == 7){
                    pout.println("ANGVELS n= "+TM.nMol+" x 3 X(0),Y(0),Z(0)...X(n-1),Y(n-1),Z(n-1)");
                    for(int j = 0; j < TM.nMol; j++){
                       pout.println(j+" "+TM.wvx[j]+" "+TM.wvy[j]+" "+TM.wvz[j]);
                    }
                }else if(i == 8){
                    pout.println("ANGACCELS n= "+TM.nMol+" x 3 X(0),Y(0),Z(0)...X(n-1),Y(n-1),Z(n-1)");
                    for(int j = 0; j < TM.nMol; j++){
                       pout.println(j+" "+TM.wax[j]+" "+TM.way[j]+" "+TM.waz[j]);
                    }
                }
            }
    }
    public void writeRestart(TheMatrix TM){
            super.openOutputWriter(file_out);
            pout = super.PW;
            for(int i = 0; i < 12; i++){
                if(i == 0){
                    pout.println("STEP "+TM.stepCount+" "+TM.restartCount);
                }else if(i == 1){
                    pout.println("MOLECULES "+TM.nMol+" "+TM.molType);
                }else if(i == 2){
                    pout.println("BOX "+(TM.regionX*TM.ro)+" "+(TM.regionY*TM.ro)+" "+(TM.regionZ*TM.ro));
                }else if(i == 3){
                    pout.println("COORDS n= "+TM.nMol+" x 3 X(0),Y(0),Z(0)...X(n-1),Y(n-1),Z(n-1)");
                    for(int j = 0; j < TM.nMol; j++){
                        pout.println(j+" "+(TM.rx[j]*TM.ro)+" "+(TM.ry[j]*TM.ro)+" "+(TM.rz[j]*TM.ro));
                    }
                }else if(i == 4){
                    pout.println("VELS n= "+TM.nMol+" x 3 X(0),Y(0),Z(0)...X(n-1),Y(n-1),Z(n-1)");
                    for(int j = 0; j < TM.nMol; j++){
                       pout.println(j+" "+TM.rvx[j]+" "+TM.rvy[j]+" "+TM.rvz[j]);
                    }
                }else if(i == 5){
                    pout.println("ACCELS n= "+TM.nMol+" x 3 X(0),Y(0),Z(0)...X(n-1),Y(n-1),Z(n-1)");
                    for(int j = 0; j < TM.nMol; j++){
                       pout.println(j+" "+TM.rax[j]+" "+TM.ray[j]+" "+TM.raz[j]);
                    }
                }else if(i == 6){
                    pout.println("ANGCOORDS n= "+TM.nMol+" x 4 q1(0),q2(0),q3(0),q(4)....q1(n-1),q2(n-1),q3(n-1),q4(n-1)");
                    for(int j = 0; j < TM.nMol; j++){
                       pout.print(j+" "+TM.rMatT[j*9]+" "+TM.rMatT[j*9+1]+" "+TM.rMatT[j*9+2]);
                       pout.print(" "+TM.rMatT[j*9+3]+" "+TM.rMatT[j*9+4]+" "+TM.rMatT[j*9+5]);
                       pout.println(" "+TM.rMatT[j*9+6]+" "+TM.rMatT[j*9+7]+" "+TM.rMatT[j*9+8]);
                    }
                }else if(i == 7){ 
                    pout.println("ANGVELS n= "+TM.nMol+" x 3 X(0),Y(0),Z(0)...X(n-1),Y(n-1),Z(n-1)");
                    for(int j = 0; j < TM.nMol; j++){
                       pout.println(j+" "+TM.wvx[j]+" "+TM.wvy[j]+" "+TM.wvz[j]);
                    }   
                }else if(i == 8){ 
                    pout.println("ANGACCELS n= "+TM.nMol+" x 3 X(0),Y(0),Z(0)...X(n-1),Y(n-1),Z(n-1)");
                    for(int j = 0; j < TM.nMol; j++){
                       pout.println(j+" "+TM.wax[j]+" "+TM.way[j]+" "+TM.waz[j]);
                    }   
                }   
            }   
    }
    public void ReadRestartHeader(TheMatrix TM){
         try{
            String line;
            String temp;
            FileReader FR;
            BufferedReader dFile;
            FR =  new FileReader(TM.restartIn);
            dFile = new BufferedReader( FR );
            int index = 0;
            int nmol = 0;
            while( ( line = dFile.readLine( ) ) != null ){
                StringTokenizer token = new StringTokenizer(line);
                StringTokenizer token2;
                temp = token.nextToken();
                if(temp.equals("STEP")){
                   TM.stepCount = Integer.parseInt(token.nextToken());
                   TM.stepLimit = TM.stepLimit + TM.stepCount;
                   TM.restartCount = Integer.parseInt(token.nextToken());
                }else if(temp.equals("BOX")){
                   TM.regionX = Double.parseDouble(token.nextToken())/TM.ro;
                   TM.regionY = Double.parseDouble(token.nextToken())/TM.ro;
                   TM.regionZ = Double.parseDouble(token.nextToken())/TM.ro;
                   TM.volume = TM.regionX*TM.regionY*TM.regionZ;
                }else if(temp.equals("MOLECULES")){
                   TM.nMol = Integer.parseInt(token.nextToken());
                   parameters.DefinePar(TM,token.nextToken());
                }else{}
             }
        }catch( IOException e ){System.err.println( e );}
    }
    public void ReadRestartFile(TheMatrix TM){
         try{
            String line;
            String temp;
            FileReader FR;
            BufferedReader dFile;
            FR =  new FileReader(TM.restartIn);
            dFile = new BufferedReader( FR );
            int index = 0;
            int nmol = 0;
            while( ( line = dFile.readLine( ) ) != null ){
                StringTokenizer token = new StringTokenizer(line);
                StringTokenizer token2;
                temp = token.nextToken();
                 if(temp.equals("COORDS")){
                    temp = token.nextToken();
                    nmol = Integer.parseInt(token.nextToken());
                    for(int i = 0; i < TM.nMol; i++){
                        line = dFile.readLine( );
                        token2 = new StringTokenizer(line);
                           temp = token2.nextToken();
                           TM.rx[i] = Double.parseDouble(token2.nextToken())/TM.ro;
                           TM.ry[i] = Double.parseDouble(token2.nextToken())/TM.ro;
                           TM.rz[i] = Double.parseDouble(token2.nextToken())/TM.ro;
                    }
                }else if(temp.equals("VELS")){
                    temp = token.nextToken();
                    nmol = Integer.parseInt(token.nextToken());
                    for(int i = 0; i < TM.nMol; i++){
                        line = dFile.readLine( );
                        token2 = new StringTokenizer(line);
                        temp = token2.nextToken();
                        TM.rvx[i] = Double.parseDouble(token2.nextToken());
                        TM.rvy[i] = Double.parseDouble(token2.nextToken());
                        TM.rvz[i] = Double.parseDouble(token2.nextToken());
                    }
                }else if(temp.equals("ACCELS")){
                    temp = token.nextToken();
                    nmol = Integer.parseInt(token.nextToken());
                    for(int i = 0; i < TM.nMol; i++){
                        line = dFile.readLine( );
                        token2 = new StringTokenizer(line);
                        temp = token2.nextToken();
                        TM.rax[i] = Double.parseDouble(token2.nextToken());
                        TM.ray[i] = Double.parseDouble(token2.nextToken());
                        TM.raz[i] = Double.parseDouble(token2.nextToken());
                    }
                }else if(temp.equals("ANGCOORDS")){
                    temp = token.nextToken();
                    nmol = Integer.parseInt(token.nextToken());
                    for(int i = 0; i < TM.nMol; i++){
                        line = dFile.readLine( );
                        token2 = new StringTokenizer(line);
                        temp = token2.nextToken();

                        TM.rMatT[i*9+0] = Double.parseDouble(token2.nextToken());
                        TM.rMatT[i*9+1] = Double.parseDouble(token2.nextToken());
                        TM.rMatT[i*9+2] = Double.parseDouble(token2.nextToken());
                        TM.rMatT[i*9+3] = Double.parseDouble(token2.nextToken());
                        TM.rMatT[i*9+4] = Double.parseDouble(token2.nextToken());
                        TM.rMatT[i*9+5] = Double.parseDouble(token2.nextToken());
                        TM.rMatT[i*9+6] = Double.parseDouble(token2.nextToken());
                        TM.rMatT[i*9+7] = Double.parseDouble(token2.nextToken());
                        TM.rMatT[i*9+8] = Double.parseDouble(token2.nextToken());
                    }
                    index = 0;
                }else if(temp.equals("ANGVELS")){
                    temp = token.nextToken();
                    nmol = Integer.parseInt(token.nextToken());
                    for(int i = 0; i < TM.nMol; i++){
                        line = dFile.readLine( );
                        token2 = new StringTokenizer(line);
                        temp = token2.nextToken();
                        TM.wvx[i] = Double.parseDouble(token2.nextToken());
                        TM.wvy[i] = Double.parseDouble(token2.nextToken());
                        TM.wvz[i] = Double.parseDouble(token2.nextToken());
                    }
                }else if(temp.equals("ANGACCELS")){
                    temp = token.nextToken();
                    nmol = Integer.parseInt(token.nextToken());
                    for(int i = 0; i < TM.nMol; i++){
                        line = dFile.readLine( );
                        token2 = new StringTokenizer(line);
                        temp = token2.nextToken();
                        TM.wax[i] = Double.parseDouble(token2.nextToken());
                        TM.way[i] = Double.parseDouble(token2.nextToken());
                        TM.waz[i] = Double.parseDouble(token2.nextToken());
                    }
                }
            }
         }catch( IOException e ){System.err.println( e );}
     }

     public void RandomGenerate(TheMatrix TM, long seed, String type){
        double[] p, tq;
        double s = 0;
        int k, k1, k2;
        p = new double[10];  tq = new double[4];
        double a1, a2, a3, fT, fR;
        double eAngx, eAngy, eAngz;
        double wvBx, wvBy, wvBz;
        double xrt, yrt, zrt;
        generator = new Random(seed);
        if(type.equals("-random")){
           for(int n = 0; n < TM.nMol; n++){ //No need to devide by TM.ro. regionXYZ already rescaled.
              TM.rx[n] = generator.nextDouble()*TM.regionX;
              TM.ry[n] = generator.nextDouble()*TM.regionY;
              TM.rz[n] = generator.nextDouble()*TM.regionZ;
           }
        }else if(type.equals("-crystal")){
           int n = 0;
           double gapX = TM.regionX/TM.initUcellX;
           double gapY = TM.regionY/TM.initUcellY;
           double gapZ = TM.regionZ/TM.initUcellZ;
           for(int nx = 0; nx < TM.initUcellX; nx++){
               for(int ny = 0; ny < TM.initUcellY; ny++){
                   for(int nz = 0; nz < TM.initUcellZ; nz++){ 
                       TM.rx[n] = nx + 0.5;  TM.rx[n] = TM.rx[n]*gapX;
                       TM.ry[n] = ny + 0.5;  TM.ry[n] = TM.ry[n]*gapY;
                       TM.rz[n] = nz + 0.5;  TM.rz[n] = TM.rz[n]*gapZ;
                       n++;
                   }
               } 
           }
        }
        for(int n = 0; n < TM.nMol; n++){
             // gen Velocities
             if(generator.nextDouble() < 0.5){TM.rvx[n] = -1*generator.nextDouble();}else{TM.rvx[n] = generator.nextDouble();}
             if(generator.nextDouble() < 0.5){TM.rvy[n] = -1*generator.nextDouble();}else{TM.rvy[n] = generator.nextDouble();}
             if(generator.nextDouble() < 0.5){TM.rvz[n] = -1*generator.nextDouble();}else{TM.rvz[n] = generator.nextDouble();}
             // gen Accelerations
             TM.rax[n] = TM.ray[n] = TM.raz[n] = 0.0;
             // gen AngCoordinates
             if(generator.nextDouble() < 0.5){xrt = -1*generator.nextDouble();}else{xrt = generator.nextDouble();}
             if(generator.nextDouble() < 0.5){yrt = -1*generator.nextDouble();}else{yrt = generator.nextDouble();}
             if(generator.nextDouble() < 0.5){zrt = -1*generator.nextDouble();}else{zrt = generator.nextDouble();}
             eAngx = Math.atan2(xrt,yrt);
             eAngy = Math.acos(zrt);
             eAngz = 2*generator.nextDouble()*Math.PI;
             a1 = 0.5 * eAngy;
             a2 = 0.5 * (eAngx - eAngz);
             a3 = 0.5 * (eAngx + eAngz);
             TM.q_u1[n] = Math.sin(a1) * Math.cos(a2);
             TM.q_u2[n] = Math.sin(a1) * Math.sin(a2);
             TM.q_u3[n] = Math.cos(a1) * Math.sin(a3);
             TM.q_u4[n] = Math.cos(a1) * Math.cos(a3);

             tq[0] = TM.q_u1[n];  tq[1] = TM.q_u2[n];  tq[2] = TM.q_u3[n];    tq[3] = TM.q_u4[n];
             for(k = 0, k2 = 0; k2 < 4; k2++){
                 for(k1 = k2; k1 < 4; k1++, k++){
                     p[k] = 2.0*tq[k1]*tq[k2];
                 }
             }
             TM.rMatT[n*9+0] = p[0] + p[9] - 1;   TM.rMatT[n*9+4] = p[4] + p[9] - 1;   TM.rMatT[n*9+8] = p[7] + p[9] - 1;
             s = 1.0;    //Transpose = 1
             TM.rMatT[n*9+1] = p[1] + s * p[8];   TM.rMatT[n*9+3] = p[1] - s * p[8];   TM.rMatT[n*9+2] = p[2] - s * p[6];
             TM.rMatT[n*9+6] = p[2] + s * p[6];   TM.rMatT[n*9+5] = p[5] + s * p[3];   TM.rMatT[n*9+7] = p[5] - s * p[3];
             //gen AngVelocities
             if(generator.nextDouble() < 0.5){TM.wvx[n] = -1*generator.nextDouble();}else{TM.wvx[n] = generator.nextDouble();}
             if(generator.nextDouble() < 0.5){TM.wvy[n] = -1*generator.nextDouble();}else{TM.wvy[n] = generator.nextDouble();}
             if(generator.nextDouble() < 0.5){TM.wvz[n] = -1*generator.nextDouble();}else{TM.wvz[n] = generator.nextDouble();}
             TM.wax[n] = TM.way[n] = TM.waz[n] = 0.0;
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
    private void command_description(){
            System.out.println("   -random     It generates a restart file from scratch. Box dimensions or side length must be entered explicitly.");
            System.out.println("               Water molecules are placed randomly in the box and it requires minimization.");
            System.out.println("               Follow the -random command by");
            System.out.println("               1) A file name for the restart file to be created in current directory.");
            System.out.println("               2) Number of water molecules.");
            System.out.println("               3) Water molecule type: tip3p, tip4p, tip5p, ST2.");
            System.out.println("               4) Leng of X box side.");
            System.out.println("               5) Leng of Y box side.");
            System.out.println("               6) Leng of Z box side.");
            System.out.println("               7) a seed for random generation of coordinates(from -2^63 to 2^63-1). This will allow reproducibility of initial configurations.");
            System.out.println("               -- All the atoms in the molecule are inside the box. This is specially important for simulations with walls.");
            System.out.println("               -- Rotational and translational Velocities are asigned random values between -1 and 1 for XYZ. Scaling to");
            System.out.println("                  the right temperatures will be done before the simulation's equilibration or production run.");
            System.out.println("               -- All randomly generated systems should be followed by minimization before production run.");
            System.out.println("   -crystal    It generates a restart file from scratch using deisred density(gm/cm^3) and number of waters only.");
            System.out.println("               unitCellX x unitCellY x unitCellZ must be equal to number of water molecules for this approach to make sense.");
            System.out.println("               Follow the -crystal command by");
            System.out.println("               1) A file name for the restart file to be created in the current directory.");
            System.out.println("               2) Number of water molecules.");
            System.out.println("               3) Molecule type: tip3p, tip4p, tip5p, ST2.");
            System.out.println("               4) Density in gm/cm^3.");
            System.out.println("               5) unitCellX.");
            System.out.println("               6) unitCellY.");
            System.out.println("               7) unitCellZ.");
            System.out.println("               8) a seed for random generation of coordinates(from -2^63 to 2^63-1). This will allow reproducibility of initial configurations.");
            System.out.println("               -- All the atoms in the molecule are inside the box. This is specially important for simulations with walls.");
            System.out.println("               -- Rotational and translational Velocities are asign random values between -1 and 1 for XYZ. Scaling to");
            System.out.println("                  the right temperatures will be done before the simulation's equilibration or production run.");
            System.out.println("               -- All randomly generated systems should be followed by minimization before production run.");
            System.out.println("   -psf        It generates a pdf file from a rst file.");
            System.out.println("               1) The restart file name.");
            System.out.println("   -crd        It generates a crd file from a rst file.");
            System.out.println("               1) The restart file name.");
            System.out.println("   -gro        It generates a gro file from a rst file.");
            System.out.println("               1) The restart file name.");
            System.out.println("   -rescale    It changes a restart file to another one with a different density. It only rescales the positions");
            System.out.println("               of the water molecules to fit in a box of a different size. Orientations, velocities and accelerations remain intact.");
            System.out.println("               Follow the -rescale command by:");
            System.out.println("               1) A path and file name of the initial restart file.");
            System.out.println("               2) A path and new file name for the rescaled restart file.");
            System.out.println("               3) A new density in gm/cm^3.");
    }
}
