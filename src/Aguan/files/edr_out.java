/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package Aguan.files;

import Aguan.TheMatrix.TheMatrix;
import java.io.*;
import java.util.*;
/**
 *
 * @author Noel Carrascal
 */
public class edr_out {
    private FileReader edr_FR, trr_FR;
    private BufferedReader edr_dataFile, trr_dataFile;
    private StringTokenizer token, token2, token3;
    public String line, line2,edr_file_name, trr_file_name;
    public String[] columnNames, columnUnits;
    public String[][] tableValues, dataType; 
    public int frameCount, index;
    public file edrOut;
    public boolean meta;
     
    public edr_out(String args1, int args3){
        try{
              String temp, temp2, temp3;
              edr_FR =  new FileReader(args1);
              edr_dataFile = new BufferedReader( edr_FR );
              index = 0;   
              int numTokens, numTokens2;
              meta = true;
              columnNames = new String[50];
              columnUnits = new String[50];    
              tableValues = new String[args3][50];
              dataType = new String[args3][50];
              while( ( line = edr_dataFile.readLine( ) ) != null ){
//              System.out.println("-"+line);
                   token = new StringTokenizer(line);
                   numTokens = token.countTokens();
                   if(numTokens > 0){
                      temp = token.nextToken();
                      if(meta){
                          if(temp.equals("gmxdump:")){ 
                             edr_file_name = token.nextToken(); 
                          }else if(temp.equals("energy")){
                             /////////////////////////////////////////////////////////////////////////////////
                             while( ( line2 = edr_dataFile.readLine( ) ) != null && meta){
                             //System.out.println("--"+line2);
                                 token2 = new StringTokenizer(line2);
                                 numTokens2 = token2.countTokens();
                                 if(numTokens2 > 0){
                                    temp2 = token2.nextToken();
                                    if(temp2.equals("time:")){
                                       meta = false;
                                       // Here, I could read the Components for the first frame.
                                       // That is unnecessary because those values are already in the edr.out file.
                                       // All I need is the Components' names and dimensions that are obtained
                                       // in the 'else' section of this 'if' statement.
                                    }else{ 
                                       index = Integer.parseInt(temp2);
                                       for(int i = 1; i < numTokens2; i++){
                                           if(i == 1){
                                              columnNames[index] = token2.nextToken();
                                           }else if(i == (numTokens2-1)){
                                              columnUnits[index] = token2.nextToken();  
                                           }else{
                                              columnNames[index] = columnNames[index]+" "+token2.nextToken();
                                           }
                                       }
                                    }
                                 }
                             }
                             /////////////////////////////////////////////////////////////////////////////////
                          }
                      }else{
                          // If, you were to obtaine all the Component values for all frames, but the first one,
                          // it would be in this 'else' section.
                      }
                   }
              } 
           edrOut = new file("./");
           edrOut.openOutputWriter("edr.out");
           edrOut.PW.printf("gmxdump: %s\n",edr_file_name);//.substring(0,eo.edr_file_name.length()));
           edrOut.PW.printf("energy components:\n");
           for(int h = 0; h < index+1; h++){
               edrOut.PW.printf("%5d  %s",h,columnNames[h]);
               for(int g = 0; g < (25-columnNames[h].length()); g++){edrOut.PW.printf(" ");}
               edrOut.PW.printf("%s\n",columnUnits[h]);
           }
           edrOut.PW.printf("\n"); 
           }catch( IOException e ){System.err.println( e );}
     }
  
     public void displayEDR(TheMatrix TM, double time, int step){
            edrOut.PW.printf("\n                   time:   %1.5e         step:%14d\n",time,step);
            edrOut.PW.printf("                                             nsteps:             1\n");
            edrOut.PW.printf("                delta_t:   2.00000e-03    sum steps:             0\n");
            edrOut.PW.printf("               Component        Energy    Av. Energy    Sum Energy\n");
            for(int hh = 0; hh < index+1; hh++){
               for(int g = 0; g < (24-columnNames[hh].length()); g++){edrOut.PW.printf(" ");}
               if(hh == 0){edrOut.PW.printf("%s  % -1.5e\n",columnNames[hh],(TM.ep*TM.uSumVDW));} // CHanged 5 for 18 to get max accuracy.
               else if(hh == 1){edrOut.PW.printf("%s  % -1.5e\n",columnNames[hh],(TM.ep*TM.uSumEE));}
               //else if(hh == 2){edrOut.PW.printf("%s % -1.5e\n",columnNames[hh],(TM.ep*TM.uSumRF1));}
               else if(hh == 2){edrOut.PW.printf("%s  % -1.5e\n",columnNames[hh],(TM.ep*TM.uSum));}
               else if(hh == 3){edrOut.PW.printf("%s  % -1.5e\n",columnNames[hh],(TM.trzKinEnergyVal));}
               else if(hh == 4){edrOut.PW.printf("%s  % -1.5e\n",columnNames[hh],((TM.ep*TM.uSum)+TM.trzKinEnergyVal));}
               else if(hh == 5){edrOut.PW.printf("%s  % -1.5e\n",columnNames[hh],TM.translationalTemperature);}
               else if(hh == 6){edrOut.PW.printf("%s  % -1.5e\n",columnNames[hh],(337.92*TM.pressure));}
               else{edrOut.PW.printf("%s  % -1.5e\n",columnNames[hh],0.0);}
            }
     }
}
