package Aguan.files;
import Aguan.TheMatrix.TheMatrix;
import java.io.*;
import java.text.*;

public class log extends file {
    public PrintWriter lout; 
    public log(String log_file){
         super(log_file);
         super.file_type = "log";
         super.openOutputWriter(log_file);
         lout = super.PW;
    }
    public void HistogramWaters(TheMatrix TM, double timeNow){
        double distHHx, distHHy, distHHz;
        String sx, sy, sz;
        TM.q = TM.q0 = TM.q1 = TM.q2 = TM.q3 = TM.q4 = TM.q5 = TM.q10 = 0;
        TM.q21 = TM.q30 = TM.q32 = TM.q40 = TM.q41 = TM.q42 = TM.q43 = TM.q50 = 0;
        TM.q51 = TM.q52 = TM.q53 = TM.q410 = TM.q421 = TM.q430 = TM.q432 = 0;
        TM.q510 = TM.q521 = TM.q530 = TM.q532 = 0;
        double[] rx2, ry2, rz2;
        rx2 = new double[TM.sitesMol];   ry2 = new double[TM.sitesMol];
        rz2 = new double[TM.sitesMol];
        double tx, ty, tz;
        double tvx, tvy, tvz, rvx, rvy, rvz;
        tvx = tvy = tvz = rvx = rvy = rvz = 0;
        NumberFormat formatter = new DecimalFormat("0.000E0");
        for(int n = 0; n < TM.nMol; n++){
        //    lout.printf("mol %d \n",n);
        //    System.out.println(n+" "+TM.rvx[n]+" "+TM.rvy[n]+" "+TM.rvz[n]);
            tvx = tvx + TM.rvx[n];   tvy = tvy + TM.rvy[n];  tvz = tvz + TM.rvz[n];
            rvx = rvx + TM.wvx[n];   rvy = rvy + TM.wvy[n];  rvz = rvz + TM.wvz[n];
            for(int m = 0; m < TM.sitesMol; m++){
                tx = TM.rMatT[(n*9)+0]*TM.rmx[m] + TM.rMatT[(n*9)+3]*TM.rmy[m] + TM.rMatT[(n*9)+6]*TM.rmz[m];
                ty = TM.rMatT[(n*9)+1]*TM.rmx[m] + TM.rMatT[(n*9)+4]*TM.rmy[m] + TM.rMatT[(n*9)+7]*TM.rmz[m];
                tz = TM.rMatT[(n*9)+2]*TM.rmx[m] + TM.rMatT[(n*9)+5]*TM.rmy[m] + TM.rMatT[(n*9)+8]*TM.rmz[m];
                rx2[m] = tx;  ry2[m] = ty;   rz2[m] = tz;
      //          lout.printf("     rx2[%d]=%3.3f   ry2[%d]=%3.3f   rz2[%d]=%3.3f\n",m,rx2[m],m,ry2[m],m,rz2[m]);
      //          lout.printf("      rx[%d]=%3.3f    ry[%d]=%3.3f    rz[%d]=%3.3f\n",m,((TM.rx[n]+tx)*TM.ro),m,((TM.ry[n]+ty)*TM.ro),m,((TM.rz[n]+tz)*TM.ro));
            }
            if(TM.molType.equalsIgnoreCase("tip3p") || TM.molType.equalsIgnoreCase("tip4p") || TM.molType.equalsIgnoreCase("spc")){
               distHHx = (rx2[2]+rx2[3])/2;    distHHy = (ry2[2]+ry2[3])/2;   distHHz = (rz2[2]+rz2[3])/2;
            }else if(TM.molType.equalsIgnoreCase("tip5p") || TM.molType.equalsIgnoreCase("st2")){
               distHHx = (rx2[2]+rx2[4])/2;    distHHy = (ry2[2]+ry2[4])/2;   distHHz = (rz2[2]+rz2[4])/2;
            }else{
               distHHx = -1;    distHHy = -1;   distHHz = -1;
               System.out.println("ERROR: Unrecognized type of molecule inside HistogramWaters.");
            }
              TM.q = getQuadrant(distHHx,distHHy,distHHz);
      //        PM.lout.printf("     distHHx=%3.3f    distHHy=%3.3f    distHHz=%3.3f  = %d\n",distHHx,distHHy,distHHz,TM.q);
                    if(TM.q == 0){TM.q0++;
              }else if(TM.q == 1){TM.q1++;
              }else if(TM.q == 2){TM.q2++;
              }else if(TM.q == 3){TM.q3++;
              }else if(TM.q == 4){TM.q4++;
              }else if(TM.q == 5){TM.q5++;
              }else if(TM.q == 10){TM.q10++;
              }else if(TM.q == 21){TM.q21++;
              }else if(TM.q == 30){TM.q30++;
              }else if(TM.q == 32){TM.q32++;
              }else if(TM.q == 40){TM.q40++;
              }else if(TM.q == 41){TM.q41++;
              }else if(TM.q == 42){TM.q42++;
              }else if(TM.q == 43){TM.q43++;
              }else if(TM.q == 50){TM.q50++;
              }else if(TM.q == 51){TM.q51++;
              }else if(TM.q == 52){TM.q52++;
              }else if(TM.q == 53){TM.q53++;
              }else if(TM.q == 410){TM.q410++;
              }else if(TM.q == 421){TM.q421++;
              }else if(TM.q == 430){TM.q430++;
              }else if(TM.q == 432){TM.q432++;
              }else if(TM.q == 510){TM.q510++;
              }else if(TM.q == 521){TM.q521++;
              }else if(TM.q == 530){TM.q530++;
              }else if(TM.q == 532){TM.q532++;
              }else{System.out.println("Error: Invalid Quadrant");}
        }
        double t0 = TM.q0 + ((TM.q10+TM.q30+TM.q40+TM.q50)/2) + ((TM.q410+TM.q430+TM.q510+TM.q530)/3);
        double t1 = TM.q1 + ((TM.q10+TM.q21+TM.q41+TM.q51)/2) + ((TM.q410+TM.q421+TM.q510+TM.q521)/3);
        double t2 = TM.q2 + ((TM.q21+TM.q32+TM.q42+TM.q52)/2) + ((TM.q421+TM.q432+TM.q521+TM.q532)/3);
        double t3 = TM.q3 + ((TM.q30+TM.q32+TM.q43+TM.q53)/2) + ((TM.q430+TM.q432+TM.q530+TM.q532)/3);
        double t4 = TM.q4 + ((TM.q40+TM.q41+TM.q42+TM.q43)/2) + ((TM.q410+TM.q421+TM.q430+TM.q432)/3);
        double t5 = TM.q5 + ((TM.q50+TM.q51+TM.q52+TM.q53)/2) + ((TM.q510+TM.q521+TM.q530+TM.q532)/3);
        sx = formatter.format(tvx); sy = formatter.format(tvy); sz = formatter.format(tvz);
        lout.printf("%5.4f %.1f %.1f %.1f %.1f %.1f %.1f ",(timeNow*1.6),t0,t1,t2,t3,t4,t5);
        lout.printf("%s %s %s %.2f %.2f %.2f \n",sx,sy,sz,rvx,rvy,rvz); 
        //lout.printf(" %d %d %d %d %d %d",TM.q10,TM.q21,TM.q30,TM.q32,TM.q40,TM.q41);
        //lout.printf(" %d %d %d %d %d %d",TM.q42,TM.q43,TM.q50,TM.q51,TM.q52,TM.q53);
        //lout.printf(" %d %d %d %d %d %d %d %d\n",TM.q410,TM.q421,TM.q430,TM.q432,TM.q510,TM.q521,TM.q530,TM.q532);
    }
    public int getQuadrant(double x, double y, double z){
        double ax, ay, az;
        ax = Math.abs(Round((float)x,3));   ay = Math.abs(Round((float)y,3));   az = Math.abs(Round((float)z,3));
              if((x>0)&&(ax>ay)&&(ax>az)){return 0;           // Volumen
        }else if((y>0)&&(ay>ax)&&(ay>az)){return 1;           // Volumen
        }else if((x<0)&&(ax>ay)&&(ax>az)){return 2;           // Volumen
        }else if((y<0)&&(ay>ax)&&(ay>az)){return 3;           // Volumen
        }else if((z>0)&&(az>ax)&&(az>ay)){return 4;           // Volumen
        }else if((z<0)&&(az>ax)&&(az>ay)){return 5;           // Volumen
        }else if((x>0)&&(y>0)&&(ax==ay)&&(ax>az)){return 10;  // Surface
        }else if((x>0)&&(y>0)&&(ax==ay)&&(ax<az)){
              if(z>0){return 4;}else if(z<0){return 5;}else{return -1;}       // Surface
        }else if((x<0)&&(y>0)&&(ax==ay)&&(ax>az)){return 21;  // Surface
        }else if((x<0)&&(y>0)&&(ax==ay)&&(ax<az)){
              if(z>0){return 4;}else if(z<0){return 5;}else{return -1;}       // Surface
        }else if((x<0)&&(y<0)&&(ax==ay)&&(ax>az)){return 32;  // Surface
        }else if((x<0)&&(y<0)&&(ax==ay)&&(ax<az)){
              if(z>0){return 4;}else if(z<0){return 5;}else{return -1;}       // Surface
        }else if((x>0)&&(y<0)&&(ax==ay)&&(ax>az)){return 30;  // Surface
        }else if((x>0)&&(y<0)&&(ax==ay)&&(ax<az)){
              if(z>0){return 4;}else if(x<0){return 5;}else{return -1;}       // Surface
        }else if((z>0)&&(x>0)&&(az==ax)&&(az>ay)){return 40;   // Surface
        }else if((z>0)&&(x>0)&&(az==ax)&&(az<ay)){
              if(y>0){return 1;}else if(y<0){return 3;}else{return -1;}       // Surface
        }else if((z>0)&&(y>0)&&(az==ay)&&(az>ax)){return 41;   // Surface
        }else if((z>0)&&(y>0)&&(az==ay)&&(az<ax)){
              if(x>0){return 0;}else if(x<0){return 2;}else{return -1;}       // Surface
        }else if((z>0)&&(x<0)&&(az==ax)&&(az>ay)){return 42;   // Surface
        }else if((z>0)&&(x<0)&&(az==ax)&&(az<ay)){
              if(y>0){return 1;}else if(y<0){return 3;}else{return -1;}       // Surface
        }else if((z>0)&&(y<0)&&(az==ay)&&(az>ax)){return 43;   // Surface
        }else if((z>0)&&(y<0)&&(az==ay)&&(az<ax)){
              if(x>0){return 0;}else if(x<0){return 2;}else{return -1;}       // Surface
        }else if((z<0)&&(x>0)&&(az==ax)&&(az>ay)){return 50;   // Surface
        }else if((z<0)&&(x>0)&&(az==ax)&&(az<ay)){
              if(y>0){return 1;}else if(y<0){return 3;}else{return -1;}       // Surface
        }else if((z<0)&&(y>0)&&(az==ay)&&(az>ax)){return 51;   // Surface
        }else if((z<0)&&(y>0)&&(az==ay)&&(az<ax)){
              if(x>0){return 0;}else if(x<0){return 2;}else{return -1;}       // Surface
        }else if((z<0)&&(x<0)&&(az==ax)&&(az>ay)){return 52;   // Surface
        }else if((z<0)&&(x<0)&&(az==ax)&&(az<ay)){
              if(y>0){return 1;}else if(y<0){return 3;}else{return -1;}       // Surface
        }else if((z<0)&&(y<0)&&(az==ay)&&(az>ax)){return 53;   // Surface
        }else if((z<0)&&(y<0)&&(az==ay)&&(az<ax)){
              if(x>0){return 0;}else if(x<0){return 2;}else{return -1;}      // Surface
        }else if((ax==ay)&&(ax==az)&&(ay==az)){
              if((x>0)&&(y>0)&&(z>0)){ return 410;
              }else if((x>0)&&(y>0)&&(z<0)){return 510;
              }else if((x>0)&&(y<0)&&(z>0)){return 430;
              }else if((x>0)&&(y<0)&&(z<0)){return 530;
              }else if((x<0)&&(y>0)&&(z>0)){return 421;
              }else if((x<0)&&(y>0)&&(z<0)){return 521;
              }else if((x<0)&&(y<0)&&(z>0)){return 432;
              }else if((x<0)&&(y<0)&&(z<0)){return 432;
              }else{return -1;}
        }else{
            return -1;
        }
    }
    public float Round(float Rval, int Rpl){
           float p = (float)Math.pow(10,Rpl);
           Rval = Rval * p;
           float tmp = Math.round(Rval);
           return (float)tmp/p;
    }
}
