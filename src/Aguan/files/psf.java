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
    private restart RR;
    public psf(String psf_file){        
           super(psf_file);
           super.openOutputWriter(psf_file);
           ppsf = super.PW;
    }
    public void createPSFFile(TheMatrix TM){
        ppsf.println("PSF CMAP");
        ppsf.println();
        ppsf.println("       2 !NTITLE");
        ppsf.println("REMARKS PSF file with "+TM.nMol+" Molecules of "+TM.molType);
        ppsf.println("REMARKS Date:##/##/##/   ##:##:##  CREATED BY USER: Noel");
        ppsf.println();
        int atomCount = 1;
        int molCount = 1;
        int particleCount, bondCount;
        particleCount = 0;
        bondCount = 0;
        //System.out.println("nMol TM.nMolNames TM.sitesMolIdx");
        for(int i = 0; i < TM.nMol; i++){
            //System.out.println(i+" "+TM.nMolNames[i]+" "+TM.sitesMolIdx[i]);
            if(TM.nMolNames[i].equals("tip5p") || TM.nMolNames[i].equals("st2")){
                particleCount = particleCount + TM.sitesMolIdx[i];
                bondCount = bondCount + 6;
            }else if(TM.nMolNames[i].equals("tip3p") || TM.nMolNames[i].equals("tip4p")){
                particleCount = particleCount + TM.sitesMolIdx[i] - 1;
                bondCount = bondCount + 3;
            }else if(TM.nMolNames[i].equals("watLJ")){
                particleCount = particleCount + TM.sitesMolIdx[i] - 1;
            }else if(TM.nMolNames[i].equalsIgnoreCase("wt0LJ") || TM.nMolNames[i].equalsIgnoreCase("wt1LJ")  || TM.nMolNames[i].equalsIgnoreCase("wt2LJ")  || TM.nMolNames[i].equalsIgnoreCase("wt3LJ")){
                particleCount = particleCount + TM.sitesMolIdx[i] - 1;
            }else if(TM.nMolNames[i].equalsIgnoreCase("wt4LJ") || TM.nMolNames[i].equalsIgnoreCase("wt5LJ")  || TM.nMolNames[i].equalsIgnoreCase("wt6LJ")  || TM.nMolNames[i].equalsIgnoreCase("wt7LJ")){
                particleCount = particleCount + TM.sitesMolIdx[i] - 1;
            }else if(TM.nMolNames[i].equalsIgnoreCase("wt8LJ") || TM.nMolNames[i].equalsIgnoreCase("wt9LJ")  || TM.nMolNames[i].equalsIgnoreCase("wtXLJ")){
                particleCount = particleCount + TM.sitesMolIdx[i] - 1;
            }
        }
        int[] bnds1 = {1,2,1,3,2,3};
        int[] bnds2 = {1,2,1,3,2,3,1,4,1,5,4,5};
        int[] bonds = new int[bondCount*2];
        int particleCount2 = 0;
        String molName;
        bondCount = 0;
        for(int i = 0; i < TM.nMol; i++){
            if(TM.nMolNames[i].equals("tip5p") || TM.nMolNames[i].equals("st2")){
                for(int j = 0; j < bnds2.length; j++){
                    bonds[bondCount] = bnds2[j] + particleCount2;
                    bondCount++;
                }
                particleCount2 = particleCount2 + TM.sitesMolIdx[i];
            }else if(TM.nMolNames[i].equals("tip3p") || TM.nMolNames[i].equals("tip4p")){
                for(int j = 0; j < bnds1.length; j++){
                    bonds[bondCount] = bnds1[j] + particleCount2;
                    bondCount++;
                }
                particleCount2 = particleCount2 + TM.sitesMolIdx[i] - 1;
            }else if(TM.nMolNames[i].equalsIgnoreCase("wt0LJ") || TM.nMolNames[i].equalsIgnoreCase("wt1LJ")  || TM.nMolNames[i].equalsIgnoreCase("wt2LJ")  || TM.nMolNames[i].equalsIgnoreCase("wt3LJ")){
                // Do nothing
            }else if(TM.nMolNames[i].equalsIgnoreCase("wt4LJ") || TM.nMolNames[i].equalsIgnoreCase("wt5LJ")  || TM.nMolNames[i].equalsIgnoreCase("wt6LJ")  || TM.nMolNames[i].equalsIgnoreCase("wt7LJ")){
                // Do nothing
            }else if(TM.nMolNames[i].equalsIgnoreCase("wt8LJ") || TM.nMolNames[i].equalsIgnoreCase("wt9LJ")  || TM.nMolNames[i].equalsIgnoreCase("wtXLJ")){
                // Do nothing
            }
        }
        ppsf.printf("%8d !NATOM\n",particleCount);
        double charge, weight;
        charge = 0.0; weight = 0.0;
        String atomTyp, atomTyp2;
        atomTyp = ""; atomTyp2 = "";
        for(int i = 0; i < TM.nMol; i++){
           if(TM.nMolNames[i].equalsIgnoreCase("tip5p")){
              molName = "TIP5";
              for(int j = 0; j < TM.sitesMolIdx[i]; j++){
                if(j == 0){
                   charge = 0.0;      atomTyp = "OH2";    atomTyp2 = "OT";    weight = 15.9994;
                }else if(j == 1){
                   charge = 0.24357;  atomTyp = "H1";     atomTyp2 = "HT";    weight = 1.00800;
                }else if(j == 2){
                   charge = 0.24357;  atomTyp = "H2";     atomTyp2 = "HT";    weight = 1.00800;
                }else if(j == 3){
                   charge = -0.24357; atomTyp = "LP1";    atomTyp2 = "HT";    weight = 0.0;
                }else if(j == 4){
                   charge = -0.24357; atomTyp = "LP2";    atomTyp2 = "HT";    weight = 0.0;
                }
                //ppsf.printf("%8d W %4d %7s %-3s %3s %12f %13.5f %11d %9.5f %18s\n",atomCount,molCount,molName,atomTyp,atomTyp2,charge,weight,0,0.0,"-0.301140E-026");
                ppsf.printf("%8d W   % -4d %5s %-3s %3s %12f %13.5f %11d %9.5f %18s\n",atomCount,molCount,molName,atomTyp,atomTyp2,charge,weight,0,0.0,"-0.301140E-026");
                atomCount++;
              }
              molCount++;
           }else if(TM.nMolNames[i].equalsIgnoreCase("st2")){
              molName = "STP2";
              for(int j = 0; j < TM.sitesMolIdx[i]; j++){
                if(j == 0){
                   charge = 0.0;      atomTyp = "OH2";    atomTyp2 = "OT";    weight = 15.9994;
                }else if(j == 1){
                   charge = 0.2410;   atomTyp = "H1";     atomTyp2 = "HT";    weight = 1.00800;
                }else if(j == 2){
                   charge = 0.2410;   atomTyp = "H2";     atomTyp2 = "HT";    weight = 1.00800;
                }else if(j == 3){
                   charge = -0.2410;  atomTyp = "LP1";    atomTyp2 = "HT";    weight = 0.0;
                }else if(j == 4){
                   charge = -0.2410;  atomTyp = "LP2";    atomTyp2 = "HT";    weight = 0.0;
                }
                //ppsf.printf("%8d W %4d %7s %-3s %3s %12f %13.5f %11d %9.5f %18s\n",atomCount,molCount,molName,atomTyp,atomTyp2,charge,weight,0,0.0,"-0.301140E-026");
                ppsf.printf("%8d W   % -4d %5s %-3s %3s %12f %13.5f %11d %9.5f %18s\n",atomCount,molCount,molName,atomTyp,atomTyp2,charge,weight,0,0.0,"-0.301140E-026");
                atomCount++;
              }
              molCount++;
           }else if(TM.nMolNames[i].equalsIgnoreCase("tip3p")){
              molName = "TIP3";
              for(int j = 0; j < (TM.sitesMolIdx[i]-1); j++){
                if(j == 0){
                   charge = -0.8340;  atomTyp = "OH2";    atomTyp2 = "OT";    weight = 15.9994;
                }else if(j == 1){
                   charge = 0.4170;   atomTyp = "H1";     atomTyp2 = "HT";    weight = 1.00800;
                }else if(j == 2){
                   charge = 0.4170;   atomTyp = "H2";     atomTyp2 = "HT";    weight = 1.00800;
                }
                //ppsf.printf("%8d W   % -4d %5s %-3s %3s %12f %13.5f %11d %9.5f %18s\n",atomCount,molCount,molName,atomTyp,atomTyp2,charge,weight,0,0.0,"-0.301140E-026");
                ppsf.printf("%8d W   % -4d %5s %-3s %3s %12f %13.5f %11d %9.5f %18s\n",atomCount,molCount,molName,atomTyp,atomTyp2,charge,weight,0,0.0,"-0.301140E-026");
                atomCount++;
              }
              molCount++;
           }else if(TM.nMolNames[i].equalsIgnoreCase("tip4p")){
              molName = "TIP4";
              for(int j = 0; j < (TM.sitesMolIdx[i]-1); j++){
                if(j == 0){
                   charge = -1.0400;  atomTyp = "OH2";    atomTyp2 = "OT";    weight = 15.9994;
                }else if(j == 1){
                   charge = 0.5200;   atomTyp = "H1";     atomTyp2 = "HT";    weight = 1.00800;
                }else if(j == 2){
                   charge = 0.5200;   atomTyp = "H2";     atomTyp2 = "HT";    weight = 1.00800;
                }
                //ppsf.printf("%8d W %4d %7s %-3s %3s %12f %13.5f %11d %9.5f %18s\n",atomCount,molCount,molName,atomTyp,atomTyp2,charge,weight,0,0.0,"-0.301140E-026");
                ppsf.printf("%8d W   % -4d %5s %-3s %3s %12f %13.5f %11d %9.5f %18s\n",atomCount,molCount,molName,atomTyp,atomTyp2,charge,weight,0,0.0,"-0.301140E-026");
                atomCount++;
              }
              molCount++;
           }else if(TM.nMolNames[i].equalsIgnoreCase("wt0LJ") || TM.nMolNames[i].equalsIgnoreCase("wt1LJ")  || TM.nMolNames[i].equalsIgnoreCase("wt2LJ")  || TM.nMolNames[i].equalsIgnoreCase("wt3LJ")){
              molName = TM.nMolNames[i].substring(2,3)+"LJ";
              for(int j = 0; j < (TM.sitesMolIdx[i]-1); j++){
                charge = 0.0;      atomTyp = "NBL";    atomTyp2 = "NB";    weight = 18.0000;
                //ppsf.printf("%8d W %4d %7s %-3s %3s %12f %13.5f %11d %9.5f %18s\n",atomCount,molCount,molName,atomTyp,atomTyp2,charge,weight,0,0.0,"-0.301140E-026");
                ppsf.printf("%8d W   % -4d %5s %-3s %3s %12f %13.5f %11d %9.5f %18s\n",atomCount,molCount,molName,atomTyp,atomTyp2,charge,weight,0,0.0,"-0.301140E-026");
                atomCount++;
              }
              molCount++;
           }else if(TM.nMolNames[i].equalsIgnoreCase("wt4LJ") || TM.nMolNames[i].equalsIgnoreCase("wt5LJ")  || TM.nMolNames[i].equalsIgnoreCase("wt6LJ")  || TM.nMolNames[i].equalsIgnoreCase("wt7LJ")){
              molName = TM.nMolNames[i].substring(2,3)+"LJ";
              for(int j = 0; j < (TM.sitesMolIdx[i]-1); j++){
                charge = 0.0;      atomTyp = "NBL";    atomTyp2 = "NB";    weight = 18.0000;
                //ppsf.printf("%8d W %4d %7s %-3s %3s %12f %13.5f %11d %9.5f %18s\n",atomCount,molCount,molName,atomTyp,atomTyp2,charge,weight,0,0.0,"-0.301140E-026");
                ppsf.printf("%8d W   % -4d %5s %-3s %3s %12f %13.5f %11d %9.5f %18s\n",atomCount,molCount,molName,atomTyp,atomTyp2,charge,weight,0,0.0,"-0.301140E-026");
                atomCount++;
              }
              molCount++;
           }else if(TM.nMolNames[i].equalsIgnoreCase("wt8LJ") || TM.nMolNames[i].equalsIgnoreCase("wt9LJ")  || TM.nMolNames[i].equalsIgnoreCase("wtXLJ")){
              molName = TM.nMolNames[i].substring(2,3)+"LJ";
              for(int j = 0; j < (TM.sitesMolIdx[i]-1); j++){
                charge = 0.0;      atomTyp = "NBL";    atomTyp2 = "NB";    weight = 18.0000;
                //ppsf.printf("%8d W %-4d %7s %-3s %3s %12f %13.5f %11d %9.5f %18s\n",atomCount,molCount,molName,atomTyp,atomTyp2,charge,weight,0,0.0,"-0.301140E-026");
                ppsf.printf("%8d W   % -4d %5s %-3s %3s %12f %13.5f %11d %9.5f %18s\n",atomCount,molCount,molName,atomTyp,atomTyp2,charge,weight,0,0.0,"-0.301140E-026");
                atomCount++;
              }
              molCount++;
           }
        }
        ppsf.println();
        ppsf.printf("%8d !NBOND\n",(bonds.length/2));
        for(int i = 0; i < bonds.length; i++){
            if((i%8)==0 && (i != 0)){
              ppsf.printf("\n");
            }
            ppsf.printf("%8d",bonds[i]);
        }
        ppsf.printf("\n\n");
        ppsf.printf("%8d !NTHETA\n",0);
        ppsf.printf("\n");
        ppsf.printf("%8d !NPHI\n",0);
        ppsf.printf("\n");
        ppsf.printf("%8d !NIMPHI\n",0);
        ppsf.printf("\n");
        ppsf.printf("%8d !NDON\n",0);
        ppsf.printf("\n");
        ppsf.printf("%8d !NACC\n",0);
    }
}
