package Aguan.files;

import Aguan.TheMatrix.TheMatrix;
import Aguan.Methods.MD;
import Aguan.parameters;
import java.io.*;
import java.util.*;
import java.text.*;

public class crd extends coordinate {
       private FileReader crd_FR;
       private BufferedReader crd_dataFile;
       public PrintWriter pcrd;
       private TheMatrix TM;
       private MD md;
       private StringTokenizer token;
       private int entryCount;
       private String line;
       private double[] fxs, fys, fzs;
       public String[] spacer = {""," ","  ","   ","    ",
                                 "     ","      ","       ",
                                 "        ","         ","          ",
                                 "           ","            ","             ","              "};
       // CRD constructor for reading from it.
       public crd(){
          super("na.crd");
          TM = new TheMatrix();
          parameters.DefinePar(TM,"tip3p");
          md = new MD();
       }
       // CRD constructor with file name for writing a CRD file
       public crd(String crd_file){
              super(crd_file);
              super.openOutputWriter(crd_file);
              pcrd = super.PW;
       }
       public void readCRD(String path){
          try{
              int tokenCount, numbTokens, lineCount, siteCount, molCount;
              String temp;
              crd_FR =  new FileReader(path);
              crd_dataFile = new BufferedReader( crd_FR );
              boolean continueLoop = true;
              entryCount = 1;
              TM.NDIM = 3;
              TM.targetTemperature = 300;
              TM.fxs = new double[TM.nMol*TM.sitesMol[0]];
              fxs = new double[TM.nMol*TM.sitesMol[0]];
              TM.fys = new double[TM.nMol*TM.sitesMol[0]];
              fys = new double[TM.nMol*TM.sitesMol[0]];
              TM.fzs = new double[TM.nMol*TM.sitesMol[0]];
              fzs = new double[TM.nMol*TM.sitesMol[0]];
              TM.rCut = 1000/TM.ro[0];
              lineCount = 0;
              siteCount = 0;
              molCount = 0;
              while(continueLoop){
                 if( ( line = crd_dataFile.readLine( ) ) != null ){
                    token = new StringTokenizer(line," \t\n\r\f,(){}[]=:");
                    numbTokens = token.countTokens();
                    temp = token.nextToken();
                    tokenCount = 0;
                    if(temp.equals("*")){
                    
                    }else{
                        if(token.countTokens() == 0){
                           TM.allAtoms = Integer.parseInt(temp);
                           TM.nMol = TM.allAtoms/3;
                           TM.initMatrix();
                           parameters.DefineMol(TM,"tip3p");
                           lineCount = -1;
                           siteCount = -1;
                        }else{
                            //count = 0;
                            for(int i = 1; i < numbTokens; i++){
                                temp = token.nextToken();
                                if(tokenCount == 2){
                                    //System.out.print(temp+" ");
                                }else if(tokenCount == 3){
                                    TM.rxs[siteCount] = Double.parseDouble(temp)/(10*TM.ro[0]);
                                    if((lineCount%3) == 0){
                                        TM.rx[molCount] = TM.rxs[siteCount];
                                        TM.rxs[siteCount+1] = TM.rxs[siteCount];
                                        //System.out.print(lineCount+" "+TM.rx[molCount]);
                                    }
                                    //System.out.print(temp+" ");
                                }else if(tokenCount == 4){
                                    TM.rys[siteCount] = Double.parseDouble(temp)/(10*TM.ro[0]);
                                    if((lineCount%3) == 0){
                                        TM.ry[molCount] = TM.rys[siteCount];
                                        TM.rys[siteCount+1] = TM.rys[siteCount];
                                        //System.out.print(" "+TM.ry[molCount]);
                                    }
                                    //System.out.print(temp+" ");
                                }else if(tokenCount == 5){
                                    TM.rzs[siteCount] = Double.parseDouble(temp)/(10*TM.ro[0]);
                                    if((lineCount%3) == 0){
                                        TM.rz[molCount] = TM.rzs[siteCount];
                                        TM.rzs[siteCount+1] = TM.rzs[siteCount];
                                        siteCount++;
                                        //System.out.println(" "+TM.rz[molCount]);
                                        molCount++;
                                    }
                                    //System.out.println(temp+" ");
                                }
                                tokenCount++;
                            }
                        }
                    }
                 }else{
                     continueLoop = false;
                 }
                 lineCount++;
                 siteCount++;
              }
              //System.out.println("Reading CRD into the Matrix");
              md.AccumProps(TM,0);
              md.ComputeSiteForces(TM);
              int counter = 1;
              double TMfxs, TMfys, TMfzs;
              for(int k = 0; k < TM.nMol*TM.sitesMol[0]; k++){
                     TM.fxs[k] = (TM.fxs[k]*TM.ep[0]/TM.ro[0]/10)/4.184;
                     TM.fys[k] = (TM.fys[k]*TM.ep[0]/TM.ro[0]/10)/4.184;
                     TM.fzs[k] = (TM.fzs[k]*TM.ep[0]/TM.ro[0]/10)/4.184;   
              }
              for(int k = 0; k < TM.nMol*TM.sitesMol[0]; k++){
                  if((k%4) == 0){
                     TMfxs = TM.fxs[k]+TM.fxs[k+1];
                     TMfys = TM.fys[k]+TM.fys[k+1];
                     TMfzs = TM.fzs[k]+TM.fzs[k+1];
                     k++;
                     //System.out.printf("%d f %.5f %.5f %.5f\n",counter,TMfxs,TMfys,TMfzs);
                     counter++;
                  }else{
                     TMfxs = TM.fxs[k];
                     TMfys = TM.fys[k];
                     TMfzs = TM.fzs[k];
                     //System.out.printf("%d f %.5f %.5f %.5f\n",counter,TMfxs,TMfys,TMfzs);
                     counter++;
                  }
              }
              //counter = 1;
              //for(int k = 0; k < TM.nMol*TM.sitesMol[0]; k++){
              //    System.out.printf("%d x %.5f %.5f %.5f\n",counter,TM.rxs[k],TM.rys[k],TM.rzs[k]);
              //    counter++;
              //}
              TM.uSum = TM.uSum*TM.ep[0]/4.184;
              TM.uSumEE = TM.uSumEE*TM.ep[0]/4.184;
              TM.uSumVDW = TM.uSumVDW*TM.ep[0]/4.184;
              //System.out.println("TM.uSum    = "+TM.uSum);
              //System.out.println("TM.uSumEE  = "+TM.uSumEE);
              //System.out.println("TM.uSumVDW = "+TM.uSumVDW);
              ComputeSiteForcesCHARMM_Correction(TM);
              //System.out.println("--------------------------------------");
              for(int k = 0; k < TM.nMol*TM.sitesMol[0]; k++){
                  if((k%4) == 0){
                     TMfxs = TM.fxs[k]+TM.fxs[k+1];
                     TMfys = TM.fys[k]+TM.fys[k+1];
                     TMfzs = TM.fzs[k]+TM.fzs[k+1];
                     k++;
                     //System.out.printf("%d f %.5f %.5f %.5f\n",counter,TMfxs,TMfys,TMfzs);
                     counter++;
                  }else{
                     TMfxs = TM.fxs[k];
                     TMfys = TM.fys[k];
                     TMfzs = TM.fzs[k];
                     //System.out.printf("%d f %.5f %.5f %.5f\n",counter,TMfxs,TMfys,TMfzs);
                     counter++;
                  }
              }
              //System.out.println("TM.uSum    = "+TM.uSum);
              //System.out.println("TM.uSumEE  = "+TM.uSumEE);
              //System.out.println("TM.uSumVDW = "+TM.uSumVDW);
           }catch( IOException e ){System.err.println( e );}
       }
       public void createCRDFile(TheMatrix TM){
          String pattern = "##0.00000";
          DecimalFormat dimentionFormatter = new DecimalFormat(pattern);
          pcrd.println("* Stepest Descend translationally Minimized waters");
          pcrd.println("* Date:##/##/##/   ##:##:##  CREATED BY USER: Noel");
          pcrd.println("*");
          //if(TM.molType.equals("tip5p") || TM.molType.equals("st2")){
          //   int numParticlesPerMol = TM.sitesMol[0];
          //   int numSizeSpace = 5 - sizeofInt(TM.nMol*numParticlesPerMol);
          //   pcrd.println(spacer[numSizeSpace]+(TM.nMol*numParticlesPerMol));
          //}else if(TM.molType.equals("tip3p") || TM.molType.equals("tip4p")){
          //   int numParticlesPerMol = TM.nCharges[0];
          //   int numSizeSpace = 5 - sizeofInt(TM.nMol*numParticlesPerMol);
          //   pcrd.println(spacer[numSizeSpace]+(TM.nMol*numParticlesPerMol));
          //}
          int atomCount = 1;
          int molCount = 0;
          if(TM.molType.equals("tip5p") || TM.molType.equals("st2")){
              molCount = 1;
          }else if(TM.molType.equals("tip3p") || TM.molType.equals("tip4p")){
              molCount = 0;
          }
          String atomTyp = "";
          String xx, yy, zz;
          int lengthx, lengthy, lengthz, lengthMC;
          boolean hydrogenFlag = true;
          boolean lonePairFlag = true;
          for(int a = 0; a < TM.nMol*TM.sitesMol[0]; a++){
              if(TM.atomType[a] != 1){
                 if(TM.atomType[a] == 2){
                    atomTyp = "OH2";
                 }else if(TM.atomType[a] == 3){
                    if(hydrogenFlag){
                       atomTyp = "H1 ";
                       hydrogenFlag = false;
                    }else{
                       atomTyp = "H2 ";
                       hydrogenFlag = true;
                    }
                 }else if(TM.atomType[a] == 4){
                    if(lonePairFlag){
                       atomTyp = "LP1";
                       lonePairFlag = false;
                    }else{
                       atomTyp = "LP2";
                       lonePairFlag = true;
                    }
                 }
                 pcrd.print(spacer[5 - sizeofInt(atomCount)]+atomCount);

                 if(TM.molType.equals("tip3p")){
                    pcrd.print(spacer[5 - sizeofInt(molCount)]+molCount+" TIP3 "+atomTyp);
                 }else if(TM.molType.equals("tip4p")){
                    pcrd.print(spacer[5 - sizeofInt(molCount)]+molCount+" TIP4 "+atomTyp);
                 }else if(TM.molType.equals("tip5p")){
                    pcrd.print(spacer[5 - sizeofInt(molCount)]+molCount+" TP5  "+atomTyp);
                 }else if(TM.molType.equals("st2")){
                    pcrd.print(spacer[5 - sizeofInt(molCount)]+molCount+" ST2  "+atomTyp);
                 }

                 xx = dimentionFormatter.format(TM.rxs[a]*TM.ro[0]);
                 lengthx = 11 - xx.length();
                 yy = dimentionFormatter.format(TM.rys[a]*TM.ro[0]);
                 lengthy = 10 - yy.length();
                 zz = dimentionFormatter.format(TM.rzs[a]*TM.ro[0]);
                 lengthz = 10 - zz.length();
                 pcrd.print(spacer[lengthx]+xx+spacer[lengthy]+yy+spacer[lengthz]+zz+" W    ");
                 pcrd.println(molCount+spacer[7 - sizeofInt(molCount)]+"0.00000");
                 atomCount++;
                 if((TM.molType.equals("tip5p") || TM.molType.equals("st2")) && atomTyp.equals("LP2")){
                     molCount++;
                 }
              }else{
                 molCount++;
              }
          }
    }
    public void createCRDFile2(TheMatrix TM){
        double tx, ty, tz;
        int index = 0;
        pcrd.println("* Stepest Descend translationally Minimized waters");
        pcrd.println("* Date:##/##/##/   ##:##:##  CREATED BY USER: Noel");
        pcrd.println("*");
        int atomCount = 1;
        int atomIndex = 0;
        int molCount = 1;
        int particleCount;
        particleCount = 0;
        for(int i = 0; i < TM.nMol; i++){
            if(TM.nMolNames[i].equals("tip5p") || TM.nMolNames[i].equals("st2")){
                particleCount = particleCount + TM.sitesMolIdx[i];
            }else if(TM.nMolNames[i].equals("tip3p") || TM.nMolNames[i].equals("tip4p")){
                particleCount = particleCount + TM.sitesMolIdx[i] - 1;
            }else if(TM.nMolNames[i].equalsIgnoreCase("wt0LJ") || TM.nMolNames[i].equalsIgnoreCase("wt1LJ") || TM.nMolNames[i].equalsIgnoreCase("wt2LJ") ||
                     TM.nMolNames[i].equalsIgnoreCase("wt3LJ") || TM.nMolNames[i].equalsIgnoreCase("wt4LJ") || TM.nMolNames[i].equalsIgnoreCase("wt5LJ") ||
                     TM.nMolNames[i].equalsIgnoreCase("wt6LJ") || TM.nMolNames[i].equalsIgnoreCase("wt7LJ") || TM.nMolNames[i].equalsIgnoreCase("wt8LJ") ||
                     TM.nMolNames[i].equalsIgnoreCase("wt9LJ") || TM.nMolNames[i].equalsIgnoreCase("wtXLJ")){
                particleCount = particleCount + TM.sitesMolIdx[i] - 1;
            }
        }
        pcrd.printf("%5d\n",particleCount);
        double charge, x, y, z;
        charge = 0.0;  x = 0.0; y = 0.0; z = 0.0;
        String atomTyp, molTyp;
        atomTyp = molTyp = "";
        for(int i = 0; i < TM.nMol; i++){
           if(TM.nMolNames[i].equals("tip5p")){
              for(int j = 0; j < TM.sitesMolIdx[i]; j++){
                if(j == 0){
                   charge = 0.0;      atomTyp = "OH2";
                }else if(j == 1){
                   charge = 0.24357;  atomTyp = "H1";
                }else if(j == 2){
                   charge = 0.24357;  atomTyp = "H2";
                }else if(j == 3){
                   charge = -0.24357; atomTyp = "LP1";
                }else if(j == 4){
                   charge = -0.24357; atomTyp = "LP2";
                }
                x = TM.rxs[index + j]*TM.ro[0];
                y = TM.rys[index + j]*TM.ro[0];
                z = TM.rzs[index + j]*TM.ro[0];
                pcrd.printf("%5d %4d %4s %-3s % 10.5f % 9.5f % 9.5f W %4d %12.5f\n",atomCount,molCount,"TIP5",atomTyp,x,y,z,molCount,charge);
                atomCount++;
                atomIndex++;
              }
              index = index + TM.sitesMolIdx[i];
              molCount++;
           }else if(TM.nMolNames[i].equals("st2")){
              for(int j = 0; j < TM.sitesMolIdx[i]; j++){
                if(j == 0){
                   charge = 0.0;      atomTyp = "OH2";
                }else if(j == 1){
                   charge = 0.2410;   atomTyp = "H1";
                }else if(j == 2){
                   charge = 0.2410;   atomTyp = "H2";
                }else if(j == 3){
                   charge = -0.2410;  atomTyp = "LP1";
                }else if(j == 4){
                   charge = -0.2410;  atomTyp = "LP2";
                }
                x = TM.rxs[index + j]*TM.ro[0];
                y = TM.rys[index + j]*TM.ro[0];
                z = TM.rzs[index + j]*TM.ro[0];
                pcrd.printf("%5d %4d %4s %-3s % 10.5f % 9.5f % 9.5f W %4d %12.5f\n",atomCount,molCount,"STP2",atomTyp,x,y,z,molCount,charge);
                atomCount++;
                atomIndex++;
              }
              index = index + TM.sitesMolIdx[i];
              molCount++;
           }else if(TM.nMolNames[i].equals("tip3p")){
              for(int j = 1; j < TM.sitesMolIdx[i]; j++){
                if(j == 1){
                   charge = -0.8340;  atomTyp = "OH2";
                }else if(j == 2){
                   charge = 0.4170;   atomTyp = "H1";
                }else if(j == 3){
                   charge = 0.4170;   atomTyp = "H2";
                }
                // x, y, z add one to j to get coords of atoms, not LJ points
                x = TM.rxs[index + j]*TM.ro[0];
                y = TM.rys[index + j]*TM.ro[0];
                z = TM.rzs[index + j]*TM.ro[0];
                pcrd.printf("%5d %4d %4s %-3s % 10.5f % 9.5f % 9.5f W %4d %12.5f\n",atomCount,molCount,"TIP3",atomTyp,x,y,z,molCount,charge);
                atomCount++;
                atomIndex++;
              }
              index = index + TM.sitesMolIdx[i];
              molCount++;
           }else if(TM.nMolNames[i].equals("tip4p")){
              for(int j = 0; j < (TM.sitesMolIdx[i]-1); j++){
                if(j == 0){
                   charge = -1.0400;  atomTyp = "OH2";
                }else if(j == 1){
                   charge = 0.5200;   atomTyp = "H1";
                }else if(j == 2){
                   charge = 0.5200;   atomTyp = "H2";
                }
                // x, y, z add one to j to get coords of atoms, not LJ points
                x = TM.rxs[index + j + 1]*TM.ro[0];
                y = TM.rys[index + j + 1]*TM.ro[0];
                z = TM.rzs[index + j + 1]*TM.ro[0];
                pcrd.printf("%5d %4d %4s %-3s % 10.5f % 9.5f % 9.5f W %4d %12.5f\n",atomCount,molCount,"TIP4",atomTyp,x,y,z,molCount,charge);
                atomCount++;
                atomIndex++;
              }
              index = index + TM.sitesMolIdx[i];
              molCount++;
           }else if(TM.nMolNames[i].equalsIgnoreCase("wt0LJ") || TM.nMolNames[i].equalsIgnoreCase("wt1LJ") || TM.nMolNames[i].equalsIgnoreCase("wt2LJ") || 
                    TM.nMolNames[i].equalsIgnoreCase("wt3LJ") || TM.nMolNames[i].equalsIgnoreCase("wt4LJ") || TM.nMolNames[i].equalsIgnoreCase("wt5LJ") || 
                    TM.nMolNames[i].equalsIgnoreCase("wt6LJ") || TM.nMolNames[i].equalsIgnoreCase("wt7LJ") || TM.nMolNames[i].equalsIgnoreCase("wt8LJ") ||
                    TM.nMolNames[i].equalsIgnoreCase("wt9LJ") || TM.nMolNames[i].equalsIgnoreCase("wtXLJ")){
              for(int j = 1; j < TM.sitesMolIdx[i]; j++){
                charge = 0.0;      atomTyp = "NBL";
                // x, y, z add one to j to get coords of atoms, not LJ points
                x = TM.rxs[index + j]*TM.ro[0];
                y = TM.rys[index + j]*TM.ro[0];
                z = TM.rzs[index + j]*TM.ro[0];
                molTyp = TM.nMolNames[i].substring(2,TM.nMolNames[i].length());
                pcrd.printf("%5d %4d %4s %-3s % 10.5f % 9.5f % 9.5f W %4d %12.5f\n",atomCount,molCount,molTyp,atomTyp,x,y,z,molCount,charge);
                atomCount++;
                atomIndex++;
              }
              index = index + TM.sitesMolIdx[i];
              molCount++;
           }
        }
    }

    public void ComputeSiteForcesCHARMM_Correction(TheMatrix TM){
        double drx, dry, drz, shiftx, shifty, shiftz;
        double fcValEE, fcValRF, rr1, rr2, rr3, rrCut2, rrCut3, rri, rri3, uVal, uValVDW, uValEE;
        double fcValx, fcValy, fcValz, virSumTemp;
        double Rmin, EPS, A, A2, A6, A12;
        int m1, m2, j1, j2, ms1, ms2, typeSum;
        rr2 = rri = 0;
        m1 = m2 = j1 = j2 = ms1 = ms2 = typeSum = 0;
        uValVDW = uValEE = fcValEE = fcValRF = 0;
        fcValx = fcValy = fcValz = virSumTemp = 0;
        rrCut2 = TM.rCut*TM.rCut;
        rrCut3 = rrCut2*TM.rCut;
        // Because this corrects for LJ interactions between H-O, and H-H in charmm that are not existing in 
        // plain Tip3 models, the forces and the energies are not zeroed. This is to prove that I can match
        // CHARMM results and that my code is valid and reliable.
        // for(int n = 0; n < TM.nMol*TM.sitesMol[0]; n++){
        //      TM.fxs[n] = 0.0; TM.fys[n] = 0.0;  TM.fzs[n] = 0.0;
        // }
        // TM.uSum = TM.uSumVDW = TM.uSumEE = TM.uSumRF1 = TM.uSumRF2 = TM.virSum = 0;
        for(m1 = 0; m1 < TM.nMol-1; m1++){
            for(m2 = m1+1; m2 < TM.nMol; m2++){
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
                    ms1 = m1 * TM.sitesMol[0];
                    ms2 = m2 * TM.sitesMol[0];
                    for(j1 = 0; j1 < TM.sitesMol[0] ; j1++){
                        for(j2 = 0; j2 < TM.sitesMol[0] ; j2++){
                            uValVDW = uValEE = fcValEE = fcValRF = 0;
                            typeSum =  TM.typeF[j1] + TM.typeF[j2];
                            if(TM.typeF[j1] == TM.typeF[j2] || typeSum == 5 || typeSum == 10){
                                drx = TM.rxs[ms1+j1] - TM.rxs[ms2+j2];
                                dry = TM.rys[ms1+j1] - TM.rys[ms2+j2];
                                drz = TM.rzs[ms1+j1] - TM.rzs[ms2+j2];
                                drx = drx + shiftx;  dry = dry + shifty;  drz = drz + shiftz;
                                drx = drx*TM.ro[0]*10;
                                dry = dry*TM.ro[0]*10;
                                drz = drz*TM.ro[0]*10;
                                // DONE BELOW SPECIALLY FOR EACH TYPE OF INTERACTION
                                //rr2 = drx*drx + dry*dry + drz*drz;
                                //rr1 = Math.sqrt(rr2);
                                //rr3 = rr2*rr1;
                                //rri = 1.0/rr2;
                                switch(typeSum){
                                    // o-o INTERACTIONS ALREADY ACCOUNTED
                                    // case 2:
                                    // rri3 = rri*rri*rri;
                                    // uValVDW = 4 * rri3 * (rri3 - 1);
                                    // fcValEE = 48 * rri3 * (rri3 - 0.5) * rri;
                                    // break;
                                    //case 4:
                                    //    uValEE =  4 * TM.bCon/rr1;
                                    //    fcValEE = uValEE * rri;
                                    //    break;
                                    case 5:
                                        // uValEE = -2 * TM.bCon/rr1;
                                        // fcValEE = uValEE * rri;
                                        rr2 = drx*drx + dry*dry + drz*drz;
                                        rr1 = Math.sqrt(rr2);
                                        Rmin = 1.768200 + 0.224500;
                                        EPS = Math.sqrt(-0.152100*(-0.046000));
                                        A = Rmin/rr1;
                                        A2 = A*A;
                                        A6 = A2*A2*A2;
                                        A12 = A6*A6;
                                        uValVDW = EPS*(A12-(2*A6));
                                        fcValEE = 12*EPS*(A6*-A12)/rr2;
                                        break;
                                    case 6:
                                        // uValEE = TM.bCon/rr1;
                                        // fcValEE = uValEE * rri;
                                        rr2 = drx*drx + dry*dry + drz*drz;
                                        rr1 = Math.sqrt(rr2);
                                        Rmin = 0.224500 + 0.224500;
                                        EPS = Math.sqrt(-0.046000*(-0.046000));
                                        A = Rmin/rr1;
                                        A2 = A*A;
                                        A6 = A2*A2*A2;
                                        A12 = A6*A6;
                                        uValVDW = EPS*(A12-(2*A6));
                                        fcValEE = 12*EPS*(A6*-A12)/rr2;
                                        break;
                                    // Considering other models other than tip3p not necessary
                                    //case 10:
                                    //    uValEE = -1 * TM.bCon/rr1;
                                    //    fcValEE = uValEE * rri;
                                    //    break;
                                    //case 14:
                                    //    uValEE = TM.bCon/rr1;
                                    //    fcValEE = uValEE * rri;
                                    //    break;
                                }
                                fcValx = (fcValEE)*drx;  fcValy = (fcValEE)*dry;     fcValz = (fcValEE)*drz;
                                TM.fxs[ms1+j1] = TM.fxs[ms1+j1] + fcValx;
                                TM.fys[ms1+j1] = TM.fys[ms1+j1] + fcValy;
                                TM.fzs[ms1+j1] = TM.fzs[ms1+j1] + fcValz;
                                
                                TM.fxs[ms2+j2] = TM.fxs[ms2+j2] - fcValx;
                                TM.fys[ms2+j2] = TM.fys[ms2+j2] - fcValy;
                                TM.fzs[ms2+j2] = TM.fzs[ms2+j2] - fcValz;

                                TM.virSum = TM.virSum + (fcValx*drx + fcValy*dry + fcValz*drz);
                                TM.uSumVDW += uValVDW;
                                TM.uSumEE += uValEE;
                                TM.uSum += (uValVDW + uValEE);
                            }
                        }
                    }
                }
            }
        }
        //System.out.print("VDW = "+(TM.uSumVDW*TM.ep[0])+",  EE = "+(TM.uSumEE*TM.ep[0]));
    }

    public int sizeofInt(int x){
           int size = 1;
           if(x > 0 && x <10) size =1;
           else if(x > 9 && x < 100)size = 2;
           else if(x > 99 && x < 1000)size = 3;
           else if(x > 999 && x < 10000)size = 4;
           else if(x > 9999 && x < 100000)size = 5;
           else if(x > 99999 && x < 1000000)size = 6;
           return size;
   }       
}
