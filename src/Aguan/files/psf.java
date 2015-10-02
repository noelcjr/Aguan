package Aguan.files;

import Aguan.TheMatrix.TheMatrix;
import java.io.*;
import java.text.*;
/**
 *
 * @author Noel Carrascal
 */
public class psf extends pointer{
    public PrintWriter ppsf;
    public String[] spacer = {""," ","  ","   ","    ",
                              "     ","      ","       ",
                              "        ","         ","          ",
                              "           ","            ","             ","              "};
    private restart RR;
    public psf(String psf_file){        
           super(psf_file);
           super.openOutputWriter(psf_file);
           ppsf = super.PW;
//           try{
  //            psfout = new FileWriter(out_psf);
    //          ppsf = new PrintWriter(psfout, true);
      //     }catch( IOException e ){System.err.println( e );}
    }
    public void createPSFFile(TheMatrix TM){
        String pattern = "##0.0000000";
        DecimalFormat dimentionFormatter = new DecimalFormat(pattern);
        String pattern2 = "##0.00000";
        DecimalFormat dimentionFormatter2 = new DecimalFormat(pattern2);
        String pattern3 = "##0.000000";
        DecimalFormat dimentionFormatter3 = new DecimalFormat(pattern3);
        String pattern4 = "##0.0000";
        DecimalFormat dimentionFormatter4 = new DecimalFormat(pattern4);
        ppsf.println("PSF CMAP");
        ppsf.println();
        ppsf.println("       2 !NTITLE");
        ppsf.println("REMARKS PSF file with "+TM.nMol+" Molecules of "+TM.molType);
        ppsf.println("REMARKS Date:##/##/##/   ##:##:##  CREATED BY USER: Noel");
        ppsf.println();
        int atomCount = 1;
        int molCount = 0;
        if(TM.molType.equals("tip5p") || TM.molType.equals("st2")){
            ppsf.println(spacer[(8 - sizeofInt(TM.nMol*TM.sitesMol))]+(TM.nMol*TM.sitesMol)+" !NATOM");
            molCount = 1;
        }else if(TM.molType.equals("tip3p") || TM.molType.equals("tip4p")){
            ppsf.println(spacer[(8 - sizeofInt(TM.nMol*TM.sitesMol))]+(TM.nMol*TM.nCharges)+" !NATOM");
            molCount = 0;
        }
        String atomTyp = "";
        String atomTyp2 = "";
        String xx, yy, zz;
        int lengthx, lengthy, lengthz;
        boolean hydrogenFlag = true;
        boolean lonePairFlag = true;
        String molName = "";
        double weight = 0;
        double[] charge;
        charge = new double[5];
        double chge = 0.0;
        int nbonds = 0;
        if(TM.molType.equals("tip3p")){
            molName = "TIP3";
            charge[0] = -0.8340;  charge[1] = 0.4170;  charge[2] = 0.4170;
            nbonds = 2;
        }else if(TM.molType.equals("tip4p")){
            molName = "TIP4";
            charge[0] = -1.0400;  charge[1] = 0.5200;  charge[2] = 0.5200;
            nbonds = 2;
        }else if(TM.molType.equals("tip5p")){
            molName = "TP5";
            charge[0] = 0.0;  charge[1] = 0.24357;  charge[2] = 0.24357;
                              charge[3] =-0.24357;  charge[4] =-0.24357;
            nbonds = 4;
        }else if(TM.molType.equals("st2")){
            molName = "ST2";
            charge[0] = 0.0;  charge[1] = 0.2410;  charge[2] = 0.2410;
                              charge[3] =-0.2410;  charge[4] =-0.2410;
            nbonds = 4;
        }else if(TM.molType.equals("spc")){
            molName = "SPC";
        }else{
            molName = "WAT";
        }
        for(int a = 0; a < TM.nMol*TM.sitesMol; a++){
            if(TM.atomType[a] != 1){
               if(TM.atomType[a] == 2){
                  atomTyp = "OH2";
                  atomTyp2 = "OT";
                  weight = 15.9994;
                  chge = charge[0];
               }else if(TM.atomType[a] == 3){
                  if(hydrogenFlag){
                     atomTyp = "H1 ";
                     atomTyp2 = "HT";
                     weight = 1.00800;
                     chge = charge[1];
                     hydrogenFlag = false;
                  }else{
                     atomTyp = "H2 ";
                     atomTyp2 = "HT";
                     weight = 1.00800;
                     chge = charge[2];
                     hydrogenFlag = true;
                  }
               }else if(TM.atomType[a] == 4){
                  if(lonePairFlag){
                     atomTyp = "LP1";
                     atomTyp2 = "HT";
                     weight = 0;
                     chge = charge[3];
                     lonePairFlag = false;
                  }else{
                     atomTyp = "LP2";
                     atomTyp2 = "HT";
                     weight = 0;
                     chge = charge[4];
                     lonePairFlag = true;                     
                  }
               }
               ppsf.print(spacer[8 - sizeofInt(atomCount)]+atomCount+" W    ");
               ppsf.print(molCount+spacer[5 - sizeofInt(molCount)]+molName+spacer[1]+atomTyp);
               ppsf.print(spacer[5-atomTyp.length()]+atomTyp2);

               xx = dimentionFormatter3.format(chge);
               if(weight == 1.00800){
                  yy = dimentionFormatter2.format(weight);
               }else if(weight == 0){
                  yy = dimentionFormatter2.format(weight);
               }else{
                  yy = dimentionFormatter4.format(weight);
               }
               ppsf.println(spacer[13 - xx.length()]+xx+"       "+yy+"           0   0.00000     -0.301140E-026");
               atomCount++;
         //      ppsf.println();
               if((TM.molType.equals("tip5p") || TM.molType.equals("st2")) && atomTyp.equals("LP2")){
                   molCount++;
               }
            }else{
               molCount++;
            }
        }
        ppsf.println();
        int column = 0;
        if(TM.molType.equals("tip3p") || TM.molType.equals("tip4p")){
           int[] bnds = {1,2,1,3,2,3};
           ppsf.println(spacer[(7 - sizeofInt(TM.nMol*bnds.length/2))]+(TM.nMol*bnds.length/2)+" !NBOND");
           for(int a = 0; a < TM.nMol; a++){
               for(int b = 0; b < bnds.length; b++){
                   ppsf.print(spacer[(7 - sizeofInt((bnds[b]+(a*TM.nCharges))))]+(bnds[b]+(a*TM.nCharges))+" ");
                   column++;
                   if(column == 8){ppsf.println(); column = 0;}
               }         
           }
           ppsf.println();
        }else if(TM.molType.equals("tip5p") || TM.molType.equals("st2")){
           int[] bnds = {1,2,1,3,2,3,1,4,1,5,4,5};
           ppsf.println(spacer[(7 - sizeofInt(TM.nMol*bnds.length/2))]+(TM.nMol*bnds.length/2)+" !NBOND");
           for(int a = 0; a < TM.nMol; a++){
               for(int b = 0; b < bnds.length; b++){
                   ppsf.print(spacer[(7 - sizeofInt((bnds[b]+(a*TM.sitesMol))))]+(bnds[b]+(a*TM.sitesMol))+" ");
                   column++;
                   if(column == 8){ppsf.println(); column = 0;}
               }           
           }
               ppsf.println();
        }
        ppsf.println();
        ppsf.println(spacer[(8 - sizeofInt(1))]+(0)+" !NTHETA");
        ppsf.println();
        ppsf.println(spacer[(8 - sizeofInt(1))]+(0)+" !NPHI");
        ppsf.println();
        ppsf.println(spacer[(8 - sizeofInt(1))]+(0)+" !NIMPHI");
        ppsf.println();
        ppsf.println(spacer[(8 - sizeofInt(1))]+(0)+" !NDON");
        ppsf.println();
        ppsf.println(spacer[(8 - sizeofInt(1))]+(0)+" !NACC");
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
