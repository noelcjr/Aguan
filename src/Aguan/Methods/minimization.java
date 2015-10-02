package Aguan.Methods;

import Aguan.TheMatrix.TheMatrix;
import Aguan.files.restart;
import Aguan.IO.dcdManager;
import Aguan.parameters;
import java.io.*;
import java.text.*;
import Aguan.files.*;
import java.util.*;
/**
 * Time Step 0.00060241
 * @author Noel Carrascal
 */
public class minimization extends MD{
    public minimization(String[] args){
        super(args);
        if(TM.boundaryCondition.equals("periodic")){
            if(TM.stepLimit > 0){
               System.out.println("System minimized for "+TM.stepLimit+" steps.");
               minJob(TM);
            }
        }else if(TM.boundaryCondition.equals("xyzWall")){
            if(TM.stepLimit > 0){
                System.out.println("System minimized for "+TM.stepLimit+" steps.");
                minJobXYZ(TM);
            }
        }else if(TM.boundaryCondition.equals("zWall")){
            if(TM.stepLimit > 0){
               System.out.println("System minimized for "+TM.stepLimit+" steps.");
               minJobZ(TM);
            }
        }else{
           System.out.println("INPUT ERROR: value for boundary condition is not valid, only");
           System.out.println("the case sensitive values of periodic, xyzWall and zWall are allowed.");
        }
        RR.writeRestart(TM);
        AccumProps(TM,0);
        try{
        dcdM.data_out.close();
        }catch( IOException e ){System.err.println( e );}
    }
    public void minJobXYZ(TheMatrix TM){
            double tfx, tfy, tfz;
            TM.stepCount = -1*TM.stepLimit;
            for(int i = 0; i < TM.stepLimit; i++){
               dcdM.write_dcdStep((TM.nPoints*TM.nMol),TM.rxs,TM.rys,TM.rzs,TM.ro,TM.vdwPoint);
               timeNow = TM.stepCount * TM.deltaT;
               AccumProps(TM, 0); 
               ComputeSiteForcesWxyz(TM);
               ComputeWallForcesXYZ(TM);
               EvalProps(TM);
               AccumProps(TM, 1); 
               OUT.PrintSummary(TM, timeNow);
               for(int a = 0; a < TM.nMol; a++){
                   tfx = tfy = tfz = 0;
                   for(int b = 0; b < TM.sitesMol ; b++){
                       tfx = tfx + TM.fxs[(a*TM.sitesMol)+b];
                       tfy = tfy + TM.fys[(a*TM.sitesMol)+b];
                       tfz = tfz + TM.fzs[(a*TM.sitesMol)+b];
                   }   
                   if(tfx > 0){TM.rx[a] = TM.rx[a] + 0.005;
                   }else{TM.rx[a] = TM.rx[a] - 0.005;}
                   if(tfy > 0){TM.ry[a] = TM.ry[a] + 0.005;
                   }else{TM.ry[a] = TM.ry[a] - 0.005;}
                   if(tfz > 0){TM.rz[a] = TM.rz[a] + 0.005;
                   }else{TM.rz[a] = TM.rz[a] - 0.005;}
               }   
               ApplyBoundaryCond(TM);
               GenSiteCoord(TM);
               TM.stepCount++;
            }   
            dcdM.write_dcdStep((TM.nPoints*TM.nMol),TM.rxs,TM.rys,TM.rzs,TM.ro,TM.vdwPoint);
    }
    public void minJobZ(TheMatrix TM){
            double tfx, tfy, tfz;
            TM.stepCount = -1*TM.stepLimit;
            for(int i = 0; i < TM.stepLimit; i++){
               dcdM.write_dcdStep((TM.nPoints*TM.nMol),TM.rxs,TM.rys,TM.rzs,TM.ro,TM.vdwPoint);
               timeNow = TM.stepCount * TM.deltaT;
               AccumProps(TM, 0); 
               ComputeSiteForcesWz(TM);
               ComputeWallForcesZ(TM);
               EvalProps(TM);
               AccumProps(TM, 1); 
               OUT.PrintSummary(TM, timeNow);
               for(int a = 0; a < TM.nMol; a++){
                   tfx = tfy = tfz = 0;
                   for(int b = 0; b < TM.sitesMol ; b++){
                       tfx = tfx + TM.fxs[(a*TM.sitesMol)+b];
                       tfy = tfy + TM.fys[(a*TM.sitesMol)+b];
                       tfz = tfz + TM.fzs[(a*TM.sitesMol)+b];
                   }   
                   if(tfx > 0){TM.rx[a] = TM.rx[a] + 0.005;
                   }else{TM.rx[a] = TM.rx[a] - 0.005;}
                   if(tfy > 0){TM.ry[a] = TM.ry[a] + 0.005;
                   }else{TM.ry[a] = TM.ry[a] - 0.005;}
                   if(tfz > 0){TM.rz[a] = TM.rz[a] + 0.005;
                   }else{TM.rz[a] = TM.rz[a] - 0.005;}
               }   
               ApplyBoundaryCond(TM);
               GenSiteCoord(TM);
               TM.stepCount++;
            }   
            dcdM.write_dcdStep((TM.nPoints*TM.nMol),TM.rxs,TM.rys,TM.rzs,TM.ro,TM.vdwPoint);
    }
    public void minJob(TheMatrix TM){
            double tfx, tfy, tfz;
            TM.stepCount = -1*TM.stepLimit;
            for(int i = 0; i < TM.stepLimit; i++){
               dcdM.write_dcdStep((TM.nPoints*TM.nMol),TM.rxs,TM.rys,TM.rzs,TM.ro,TM.vdwPoint);
               timeNow = TM.stepCount * TM.deltaT;
               AccumProps(TM, 0);
               ComputeSiteForces(TM);
               EvalProps(TM);
               AccumProps(TM, 1);
               OUT.PrintSummary(TM, timeNow);
               for(int a = 0; a < TM.nMol; a++){
                   tfx = tfy = tfz = 0;
                   for(int b = 0; b < TM.sitesMol ; b++){
                       tfx = tfx + TM.fxs[(a*TM.sitesMol)+b];
                       tfy = tfy + TM.fys[(a*TM.sitesMol)+b];
                       tfz = tfz + TM.fzs[(a*TM.sitesMol)+b];
                   }
                   if(tfx > 0){TM.rx[a] = TM.rx[a] + 0.005;
                   }else{TM.rx[a] = TM.rx[a] - 0.005;}
                   if(tfy > 0){TM.ry[a] = TM.ry[a] + 0.005;
                   }else{TM.ry[a] = TM.ry[a] - 0.005;}
                   if(tfz > 0){TM.rz[a] = TM.rz[a] + 0.005;
                   }else{TM.rz[a] = TM.rz[a] - 0.005;}
               }
               ApplyBoundaryCond(TM);
               GenSiteCoord(TM);
               TM.stepCount++;
            }
            dcdM.write_dcdStep((TM.nPoints*TM.nMol),TM.rxs,TM.rys,TM.rzs,TM.ro,TM.vdwPoint);
    }
    public void minimization(TheMatrix TM){
        double dx, dy, dz;
        double energy0 = 0.0;
        double energy1 = 0.0;
        int counter = 0;
        for(int n = 0; n < TM.nMol; n++){
            ComputeSiteForces(TM);
            energy0 = TM.uSum;
            for(int m = 0; m < 6; m++){
                translate(TM, n, m, 0.1);
                ApplyBoundaryCond(TM);
                ComputeSiteForces(TM);
                energy1 =  TM.uSum;
                if(energy0 < energy1){
                   translate(TM, n, m,-0.1);
                }else{
                   if(m == 0){
                      TM.rx[n] = TM.rxs[n*TM.sitesMol];
                   }else if(m == 1){
                      TM.rx[n] = TM.rxs[n*TM.sitesMol];
                   }else if(m == 2){
                      TM.ry[n] = TM.rys[n*TM.sitesMol];
                   }else if(m == 3){
                      TM.ry[n] = TM.rys[n*TM.sitesMol];
                   }else if(m == 4){
                      TM.rz[n] = TM.rzs[n*TM.sitesMol];
                   }else if(m == 5){
                      TM.rz[n] = TM.rzs[n*TM.sitesMol];
                   }
                }
            }
        }
   //     ComputeSiteForces(TM);  // This routine sems redundant. Check and see if it can
                                  // be eliminated.
    }
    public void translate(TheMatrix TM, int n, int dim, double d){
        if(dim == 0){
           for(int j = 0; j < TM.sitesMol; j++){
               TM.rxs[n*TM.sitesMol + j] = TM.rxs[n*TM.sitesMol + j] + d;
           }
        }else if(dim == 1){
           for(int j = 0; j < TM.sitesMol; j++){
               TM.rxs[n*TM.sitesMol + j] = TM.rxs[n*TM.sitesMol + j] - d;
           }
        }else if(dim == 2){
            for(int j = 0; j < TM.sitesMol; j++){
               TM.rys[n*TM.sitesMol + j] = TM.rys[n*TM.sitesMol + j] + d;
            }
        }else if(dim == 3){
            for(int j = 0; j < TM.sitesMol; j++){
               TM.rys[n*TM.sitesMol + j] = TM.rys[n*TM.sitesMol + j] - d;
            }
        }else if(dim == 4){
            for(int j = 0; j < TM.sitesMol; j++){
               TM.rzs[n*TM.sitesMol + j] = TM.rzs[n*TM.sitesMol + j] + d;
            }
        }else if(dim == 5){
            for(int j = 0; j < TM.sitesMol; j++){
               TM.rzs[n*TM.sitesMol + j] = TM.rzs[n*TM.sitesMol + j] - d;
            }
        }
    }
}
