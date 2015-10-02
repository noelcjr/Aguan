package Aguan.files;

import Aguan.TheMatrix.TheMatrix;
import java.io.*;
import java.text.*;

public class crd extends coordinate {
       public PrintWriter pcrd;
       public String[] spacer = {""," ","  ","   ","    ",
                              "     ","      ","       ",
                              "        ","         ","          ",
                              "           ","            ","             ","              "};
       public crd(String crd_file){
              super(crd_file);
              super.openOutputWriter(crd_file);
              pcrd = super.PW;
       }
       public void createCRDFile(TheMatrix TM){
          String pattern = "##0.00000";
          DecimalFormat dimentionFormatter = new DecimalFormat(pattern);
          pcrd.println("* Stepest Descend translationally Minimized waters");
          pcrd.println("* Date:##/##/##/   ##:##:##  CREATED BY USER: Noel");
          pcrd.println("*");
          if(TM.molType.equals("tip5p") || TM.molType.equals("st2")){
             int numParticlesPerMol = TM.sitesMol;
             int numSizeSpace = 5 - sizeofInt(TM.nMol*numParticlesPerMol);
             pcrd.println(spacer[numSizeSpace]+(TM.nMol*numParticlesPerMol));
          }else if(TM.molType.equals("tip3p") || TM.molType.equals("tip4p")){
             int numParticlesPerMol = TM.nCharges;
             int numSizeSpace = 5 - sizeofInt(TM.nMol*numParticlesPerMol);
             pcrd.println(spacer[numSizeSpace]+(TM.nMol*numParticlesPerMol));
          }

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
          for(int a = 0; a < TM.nMol*TM.sitesMol; a++){
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

                 xx = dimentionFormatter.format(TM.rxs[a]*TM.ro);
                 lengthx = 11 - xx.length();
                 yy = dimentionFormatter.format(TM.rys[a]*TM.ro);
                 lengthy = 10 - yy.length();
                 zz = dimentionFormatter.format(TM.rzs[a]*TM.ro);
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
