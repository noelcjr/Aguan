/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package Aguan.files;

import Aguan.TheMatrix.TheMatrix;
import java.io.*;
import java.util.*;
import Aguan.Methods.MD;
import Aguan.parameters;
import Aguan.files.edr_out;

/**
 *
 * @author Noel Carrascal
 */
public class trr_out {
    private FileReader trr_FR;
    private BufferedReader trr_dataFile;
    private String line, trr_file_name, temp, first_token, second_token;
    private StringTokenizer token, token2;
    private TheMatrix TM;
    private MD md;
    public String[] columnNames, columnUnits; 
    public int[] frame;
    public int frameCount, entryCount;
    public file edr_out;
    public file trr_out;
    private double[] fxs, fys, fzs; 
    private edr_out eo;

    public int natoms, step, index, lamda;
    public double time;
    public double[] fxs2, fys2, fzs2;
    public double[][] box;
    public trr_out(String[] args){
        try{
           eo = new edr_out(args[1],Integer.parseInt(args[3]));  
           trr_FR =  new FileReader(args[2]);                    
           trr_dataFile = new BufferedReader( trr_FR );         
           trr_out = new file("./");
           trr_out.openOutputWriter("trr.out");
           TM = new TheMatrix();
           md = new MD();
           parameters.DefinePar(TM,"tip3p");
           TM.NDIM = 3;
           TM.allAtoms = Integer.parseInt(args[4]);
           TM.nMol = TM.allAtoms/3;
           //  System.out.println("TM.nMol="+TM.nMol);
           //  System.out.println("TM.sitesMol[0]="+TM.sitesMol[0]);
           //  System.out.println("TM.bCon="+TM.bCon);
           //  System.out.println("TM.ro[0]="+TM.ro[0]);
           //  System.out.println("TM.ep[0]="+TM.ep[0]);
           //  System.out.println("TM.electricField="+TM.electricField);
           TM.fxs = new double[TM.nMol*TM.sitesMol[0]];
           fxs = new double[TM.nMol*TM.sitesMol[0]];
           TM.fys = new double[TM.nMol*TM.sitesMol[0]];
           fys = new double[TM.nMol*TM.sitesMol[0]];
           TM.fzs = new double[TM.nMol*TM.sitesMol[0]];
           fzs = new double[TM.nMol*TM.sitesMol[0]];
           TM.targetTemperature = Double.parseDouble(args[5]);
           TM.rCut = Double.parseDouble(args[6])/TM.ro[0];
           // System.out.println("TM.rCut="+TM.rCut+" "+args[6]+"x10 A");
           TM.tau = Double.parseDouble(args[7]);
           // System.out.println("TM.tau="+TM.tau);
           TM.initMatrix();
           parameters.DefineMol(TM,"tip3p");
           int frame = 0, index = 0;   int numTokens, numTokens2;
           box = new double[3][3];       
           boolean continueLoop = true;
           frameCount = 0;
           entryCount = 1;
           System.out.println("count frame atom dimesion GromacsF AguanF");
           while(continueLoop){
                     if( ( line = trr_dataFile.readLine( ) ) != null ){
                       token = new StringTokenizer(line," \t\n\r\f,(){}[]=:");
                       numTokens = token.countTokens();
                       first_token = token.nextToken();
                       second_token = token.nextToken();
                       if(second_token.equals("frame")){
                          trr_file_name = first_token;
                          frame = Integer.parseInt(token.nextToken());
                          //System.out.println(trr_file_name+"-->"+frame);
                       }else if(first_token.equals("natoms")){
                             natoms = Integer.parseInt(second_token);
                             if(token.nextToken().equals("step")){
                                step = Integer.parseInt(token.nextToken());
                             }
                             if(token.nextToken().equals("time")){
                                time = Double.parseDouble(token.nextToken());
                             }
                             if(token.nextToken().equals("lamda")){
                                lamda = Integer.parseInt(token.nextToken());
                             }
                       }else if(first_token.equals("box")){
                             if(second_token.equals("0")){
                                box[0][0] = Double.parseDouble(token.nextToken());
                                box[0][1] = Double.parseDouble(token.nextToken());
                                box[0][2] = Double.parseDouble(token.nextToken());
                                TM.regionX = box[0][0]/TM.ro[0];
                             }else if(second_token.equals("1")){
                                box[1][0] = Double.parseDouble(token.nextToken());
                                box[1][1] = Double.parseDouble(token.nextToken());
                                box[1][2] = Double.parseDouble(token.nextToken());
                                TM.regionY = box[1][1]/TM.ro[0];
                             }else if(second_token.equals("2")){
                                box[2][0] = Double.parseDouble(token.nextToken());
                                box[2][1] = Double.parseDouble(token.nextToken());
                                box[2][2] = Double.parseDouble(token.nextToken());
                                TM.regionZ = box[2][2]/TM.ro[0];
                             }
                             TM.volume = TM.regionX*TM.regionY*TM.regionZ;
                       }else if(first_token.equals("x")){
                             int index2 = 0;  int counter = 0;
                             for(int i = 0; i < natoms; i++){
                                 if( ( line = trr_dataFile.readLine( ) ) != null ){
                                   token2 = new StringTokenizer(line," \t\n\r\f,(){}[]=:");
                                   temp = token2.nextToken();
                                   index = Integer.parseInt(token2.nextToken());
                                      TM.rxs[counter] = Double.parseDouble(token2.nextToken())/TM.ro[0];
                                      TM.rys[counter] = Double.parseDouble(token2.nextToken())/TM.ro[0];
                                      TM.rzs[counter] = Double.parseDouble(token2.nextToken())/TM.ro[0];
                                      //System.out.printf("      x[%5d]={% -1.5e, % -1.5e, % -1.5e}\n",counter,(TM.rxs[counter]*TM.ro[0]),(TM.rys[counter]*TM.ro[0]),(TM.rzs[counter]*TM.ro[0]));
                                      if((index%3) == 0){
                                        TM.rx[index2] = TM.rxs[counter];
                                        TM.ry[index2] = TM.rys[counter];
                                        TM.rz[index2] = TM.rzs[counter];
                                        index2++;
                                        counter++;
                                        TM.rxs[counter] = TM.rxs[counter-1];
                                        TM.rys[counter] = TM.rys[counter-1];
                                        TM.rzs[counter] = TM.rzs[counter-1];
                                        //System.out.printf("      x[%5d]={% -1.5e, % -1.5e, % -1.5e}\n",counter,(TM.rxs[counter]*TM.ro[0]),(TM.rys[counter]*TM.ro[0]),(TM.rzs[counter]*TM.ro[0]));
                                      }
                                      counter++;
                                 }
                             } 
                             //System.out.println("index2="+index2);
                       }else if(first_token.equals("v")){
                             // Reads velocities and Calculates temperatures at the same time
                             TM.vvSum = 0.0;
                             for(int j = 0; j < natoms; j++){
                                 if( ( line = trr_dataFile.readLine( ) ) != null ){
                                   token2 = new StringTokenizer(line," \t\n\r\f,(){}[]=:");
                                   temp = token2.nextToken();
                                   index = Integer.parseInt(token2.nextToken());
                                   TM.vxs[index] = Double.parseDouble(token2.nextToken());
                                   TM.vys[index] = Double.parseDouble(token2.nextToken());
                                   TM.vzs[index] = Double.parseDouble(token2.nextToken());
                                    //System.out.printf("      x[%5d]={% -1.5e, % -1.5e, % -1.5e}\n",counter,TM.vxs[index],TM.vys[index],TM.vzs[index]);
                                 }
                             }
                             for(int n = 0; n <= index; n++){
                                  if(n == 0 || n ==3){
                                     TM.vvSum += 16*(TM.vxs[n]*TM.vxs[n] + TM.vys[n]*TM.vys[n] + TM.vzs[n]*TM.vzs[n]);
                                  }else{
                                     TM.vvSum += TM.vxs[n]*TM.vxs[n] + TM.vys[n]*TM.vys[n] + TM.vzs[n]*TM.vzs[n];
                                  }
                             }
                             TM.trzKinEnergyVal = 0.5 * TM.vvSum;
                             TM.translationalTemperature = TM.trzKinEnergyVal/(6*0.00831451); // T = 2 KE/Nkb
                       }else if(first_token.equals("f")){
                             for(int k = 0; k < natoms; k++){
                                 if( ( line = trr_dataFile.readLine( ) ) != null ){
                                   token2 = new StringTokenizer(line," \t\n\r\f,(){}[]=:");
                                   temp = token2.nextToken();
                                   index = Integer.parseInt(token2.nextToken());
                                   fxs[index] = Double.parseDouble(token2.nextToken());
                                   fys[index] = Double.parseDouble(token2.nextToken());
                                   fzs[index] = Double.parseDouble(token2.nextToken());
                                 }
                             }
                             md.ComputeSiteForces(TM);
                             int counter = 1;
                             for(int k = 0; k < natoms; k++){
                                 System.out.printf("%7d %d %d x % -1.5e % -1.5e\n",entryCount,frameCount,counter,TM.fxs[k],fxs[k]);
                                 entryCount++;
                                 System.out.printf("%7d %d %d y % -1.5e % -1.5e\n",entryCount,frameCount,counter,TM.fys[k],fys[k]);
                                 entryCount++;
                                 System.out.printf("%7d %d %d z % -1.5e % -1.5e\n",entryCount,frameCount,counter,TM.fzs[k],fzs[k]);
                                 entryCount++;
                                 counter++;
                             }
                             TM.pressure = (2/(TM.volume*3))*(0.5*TM.vvSum + 0.5*TM.virSum);
                             displayTRR(frame); 
                             eo.displayEDR(TM,time,step);
                             frameCount++;
                             //System.out.println("%%%%%%%%%%%%%%%%%%%%%%%%");         
                       }
                     }else{
                          continueLoop = false;
                     }
           }
        }catch( IOException e ){System.err.println( e );}
     }

    public void displayTRR(int frame){
            int c, d;
            trr_out.PW.printf("%s frame %d:\n",trr_file_name,frame);
            trr_out.PW.printf("   natoms=%10d  step=%10d  time=%1.7e  lambda=%10d\n",natoms,step,time,lamda);
            trr_out.PW.printf("   box (3x3):\n");
            trr_out.PW.printf("      box[    0]={ %1.5e,  %1.5e,  %1.5e}\n",(box[0][0]),(box[0][1]),(box[0][2]));
            trr_out.PW.printf("      box[    1]={ %1.5e,  %1.5e,  %1.5e}\n",(box[1][0]),(box[1][1]),(box[1][2]));
            trr_out.PW.printf("      box[    2]={ %1.5e,  %1.5e,  %1.5e}\n",(box[2][0]),(box[2][1]),(box[2][2]));
            trr_out.PW.println("   x ("+natoms+"x"+"3):");
            int counter = 0;
            for(c = 0; c < natoms+TM.nMol; c++){
                if((c%TM.sitesMol[0])==0){
                }else{
                     trr_out.PW.printf("      x[%5d]={% -1.5e, % -1.5e, % -1.5e}\n",counter,(TM.ro[0]*TM.rxs[c]),(TM.ro[0]*TM.rys[c]),(TM.ro[0]*TM.rzs[c]));
                     counter++;
                }
            }
            trr_out.PW.println("   v ("+natoms+"x"+"3):");
            for(c = 0; c < natoms; c++){
                trr_out.PW.printf("      v[%5d]={% -1.5e, % -1.5e, % -1.5e}\n",c,TM.vxs[c],TM.vys[c],TM.vzs[c]);
            }
            trr_out.PW.println("   f ("+natoms+"x"+"3):");
            double tempX, tempY, tempZ;
            int index = 0; int index2 = 0;
            for(c = 0; c < TM.nMol; c++){
                tempX = tempY = tempZ = 0;
                index2 = c*TM.sitesMol[0];
                for(d = 0; d < TM.sitesMol[0]; d++){
                    if(d == 0){
                       tempX = TM.fxs[index2+d];   tempY = TM.fys[index2+d];    tempZ = TM.fzs[index2+d];
                    }else if(d==1){
                       tempX = tempX + TM.fxs[index2+d];   tempY = tempY + TM.fys[index2+d];    tempZ = tempZ + TM.fzs[index2+d];
                       trr_out.PW.printf("      f[%5d]={% -1.5e, % -1.5e, % -1.5e}\n",index,(TM.ep[0]*tempX/TM.ro[0]),(TM.ep[0]*tempY/TM.ro[0]),(TM.ep[0]*tempZ/TM.ro[0]));
                       index++;
                    }else{
                       trr_out.PW.printf("      f[%5d]={% -1.5e, % -1.5e, % -1.5e} c= %d , d= %d\n",index,(TM.ep[0]*TM.fxs[index2+d]/TM.ro[0]),(TM.ep[0]*TM.fys[index2+d]/TM.ro[0]),(TM.ep[0]*TM.fzs[index2+d]/TM.ro[0]),c,d);
                       index++;
                    }
                }
            }
     }    
}
