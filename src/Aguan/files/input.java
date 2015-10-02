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
 * @author Noel
 */
public class input {
private FileReader FR;
private BufferedReader dataFile;
private String line;

public input(String input) {
        try{
           FR =  new FileReader(input);
           dataFile = new BufferedReader( FR );
        }catch( IOException e ){System.err.println( e );}
}
// This could be improved to check for all missing input entries.
// This would avoid defaults that are set without user concent.
// Even if those defaults may seem like reasonable.
public void GetNameList(TheMatrix TM){
        try{
           String temp;
           while( ( line = dataFile.readLine( ) ) != null ){
               StringTokenizer token = new StringTokenizer(line);
               temp = token.nextToken();
                     if(temp.equalsIgnoreCase("deltaT")){                TM.deltaT = Double.parseDouble(token.nextToken());
               }else if(temp.equalsIgnoreCase("boundaryCondition")){     TM.boundaryCondition = token.nextToken();
               }else if(temp.equalsIgnoreCase("rCut")){                  TM.rCut = Double.parseDouble(token.nextToken());
               }else if(temp.equalsIgnoreCase("stepAvg")){               TM.stepAvg = Integer.parseInt(token.nextToken());
               }else if(temp.equalsIgnoreCase("trajOut")){               TM.trajOut = Integer.parseInt(token.nextToken());
               }else if(temp.equalsIgnoreCase("temperature")){           TM.targetTemperature = Double.parseDouble(token.nextToken());
               }else if(temp.equalsIgnoreCase("tau")){                   TM.tau = Double.parseDouble(token.nextToken());
               }else if(temp.equalsIgnoreCase("electricField")){         TM.electricField = Double.parseDouble(token.nextToken());
               }else if(temp.equalsIgnoreCase("TempRateChange")){        TM.TempRateChange = Double.parseDouble(token.nextToken());
               }else if(temp.equalsIgnoreCase("ThermostatType")){        TM.ThermostatType = Integer.parseInt(token.nextToken());
               }else if(temp.equalsIgnoreCase("wallType")){              TM.wallType = Integer.parseInt(token.nextToken());
               }
           }
        }catch( IOException e ){System.err.println( e );}
}

public void GetNameList_pr_08_4(TheMatrix TM){
    try{
        String temp;
        while( ( line = dataFile.readLine( ) ) != null ){
            StringTokenizer token = new StringTokenizer(line);
            temp = token.nextToken();
                  if(temp.equalsIgnoreCase("deltaT")){    TM.deltaT = Double.parseDouble(token.nextToken());
            }else if(temp.equalsIgnoreCase("density")){   TM.density = Double.parseDouble(token.nextToken());
            }else if(temp.equalsIgnoreCase("initUcell")){ TM.initUcellX = Integer.parseInt(token.nextToken());
                                                          TM.initUcellY = Integer.parseInt(token.nextToken());
                                                          TM.initUcellZ = Integer.parseInt(token.nextToken());
            }else if(temp.equalsIgnoreCase("rCut")){      TM.rCut = Double.parseDouble(token.nextToken());
            }else if(temp.equalsIgnoreCase("stepAdjustTemp")){TM.stepAdjustTemp = Integer.parseInt(token.nextToken()); 
            }else if(temp.equalsIgnoreCase("stepAvg")){    TM.stepAvg = Integer.parseInt(token.nextToken());
            }else if(temp.equalsIgnoreCase("stepEquil")){  TM.stepEquil = Integer.parseInt(token.nextToken());
            }else if(temp.equalsIgnoreCase("stepLimit")){  TM.stepLimit = Integer.parseInt(token.nextToken());
            }else if(temp.equalsIgnoreCase("temperature")){TM.targetTemperature = Double.parseDouble(token.nextToken());
            }
       }
    }catch( IOException e ){System.err.println( e );}
}
}
