package Aguan.ntraj;

import java.io.*;
import java.util.*;
import Aguan.files.*;
public class diffusion {
    private PrintWriter doutC, doutV, doutF;
    private FileWriter outcrd, outvel, outfrc;
    private FileReader FR;
    private BufferedReader dFile;
    private DataInputStream dis;
    private String prefix;
    private int frame1, natom1, step1;
    private float time1, lamda1;
    private int frame2, natom2, step2;
    private float time2, lamda2;
    private float[][] box1, box2, crd1, crd2, vel1, vel2, frc1, frc2;
    private String[] frameData1, frameData2;
    private float sumX, sumY, sumZ;
    private float msd, msdX, msdY, msdZ;
    private dcd DCD;
    public diffusion(String filePath, String fileType){
       try{
           prefix = filePath.substring(0,(filePath.length()-4));
           System.out.println("Reading file trj_typ: "+prefix);
           if(fileType.equals("dcd")){
              DCD = new dcd(filePath);
              DCD.openInputStream();
              System.out.printf("DCD file %s \n",DCD.read_dcd());
              DCD.est_file_length();
              DCD.print_header();
              DCD.gen_DCD_body_pointers();
           }else if(fileType.equals("trr")){
               FR =  new FileReader(filePath);
               dFile = new BufferedReader( FR );
           }
           outcrd = new FileWriter(prefix+"_crd.xvg");
           doutC = new PrintWriter(outcrd, true);
           outvel = new FileWriter(prefix+"_vel.xvg");
           doutV = new PrintWriter(outvel, true);
           outfrc = new FileWriter(prefix+"_frc.xvg");
           doutF = new PrintWriter(outfrc, true);
           frameData1 = new String[5];   frameData2 = new String[5];
           sumX = sumY = sumZ = 0;
        }catch( IOException e ){System.err.println( e );} 
    }
    public void read_dcd_frame(int n, int nset){
        int i, j;
        for(i = 0; i < 3; i++){
            System.out.println(DCD.read_int());
            for(j = 0; j < n; j++){
                System.out.println("["+i+","+j+"]="+DCD.read_float());
            }
            System.out.println(DCD.read_int());
        }
    }
    public void read_dcd_file(){
         String line;
         int lineCount, frameCount, natoms, i, NSET, N;
         float time;
         lineCount = frameCount = 0;
         N = DCD.get_N();     NSET = DCD.get_NSET();
    }
    public void read_trr_file(){
         try{
            String line;
            int lineCount, frameCount, natoms, i;
            float time;
            lineCount = frameCount = 0;
            while( ( line = dFile.readLine( ) ) != null ){
                 if(frameCount == 0){
                    read_header(line, frameData1);
                    natoms = Integer.parseInt(frameData1[1]);
                    box1 = new float[3][3];    crd1 = new float[natoms][3];  vel1 = new float[natoms][3];  frc1 = new float[natoms][3];
                    read_frame(box1, crd1, vel1, frc1);
                    frameCount++;
                 }else if(frameCount == 1){
                    read_header(line, frameData2);
                    natoms = Integer.parseInt(frameData2[1]);
                    time = Float.parseFloat(frameData2[3]);
                    box2 = new float[3][3];    crd2 = new float[natoms][3];  vel2 = new float[natoms][3];  frc2 = new float[natoms][3];
                    read_frame(box2, crd2, vel2, frc2);
                    write_xvg_mean_sd(time, natoms, crd1, crd2, box2, "crd");
                    write_MSD(time, natoms, crd1, crd2, box2);
                    write_xvg_mean_sd(time, natoms, vel1, vel2, box2, "vel");
                    write_xvg_mean_sd(time, natoms, frc1, frc2, box2, "frc");
         //           System.out.println(frameCount+" "+time+" crd1, crd2: "+frameData1[1]+" "+frameData2[1]+" "+frameData1[2]+" "+frameData2[2]+" "+frameData1[3]+" "+frameData2[3]);
                    frameCount++;
                 }else{
                    if((frameCount%2) == 0){
                      read_header(line, frameData1);
                      natoms = Integer.parseInt(frameData1[1]);
                      time = Float.parseFloat(frameData1[3]);
                      box1 = new float[3][3];    crd1 = new float[natoms][3];  vel1 = new float[natoms][3];  frc1 = new float[natoms][3];
                      read_frame(box1,  crd1, vel1, frc1);
                      write_xvg_mean_sd(time, natoms, crd2, crd1, box1, "crd");
                      write_xvg_mean_sd(time, natoms, vel2, vel1, box1, "vel");
                      write_xvg_mean_sd(time, natoms, frc2, frc1, box1, "frc");
         //             System.out.println(frameCount+" "+time+" crd2, crd1: "+frameData2[1]+" "+frameData1[1]+" "+frameData2[2]+" "+frameData1[2]+" "+frameData2[3]+" "+frameData1[3]);
                      frameCount++;
                    }else{
                      read_header(line, frameData2);
                      natoms = Integer.parseInt(frameData2[1]);
                      time = Float.parseFloat(frameData2[3]);
                      box2 = new float[3][3];    crd2 = new float[natoms][3];  vel2 = new float[natoms][3];  frc2 = new float[natoms][3];
                      read_frame(box2, crd2, vel2, frc2);
                      write_xvg_mean_sd(time, natoms, crd1, crd2, box2, "crd");
                      write_xvg_mean_sd(time, natoms, vel1, vel2, box2, "vel");
                      write_xvg_mean_sd(time, natoms, frc1, frc2, box2, "frc");
         //             System.out.println(frameCount+" "+time+" crd1, crd2: "+frameData1[1]+" "+frameData2[1]+" "+frameData1[2]+" "+frameData2[2]+" "+frameData1[3]+" "+frameData2[3]);
                      frameCount++;
                    }
                 }
            }
            
         }catch( IOException e ){System.err.println( e );}
    }
    private void write_MSD(float fc, int r, float[][] arr1, float[][] arr2, float[][] box){
            int nOxygen = r/3;
            float[][] del = new float[nOxygen][3];
            int i, j;
            for(i = 0; i < nOxygen; i++){
                del[i][0] = arr2[i*3][0] - arr1[i*3][0];    del[i][1] = arr2[i*3][1] - arr1[i*3][1];     del[i][2] = arr2[i*3][2] - arr1[i*3][2];
                for(j = 0; j < 3; j++){
                   if((del[i][j] < 0) && (del[i][j] < (-1*box[j][j]/2))){
                       del[i][j] = del[i][j] + box[j][j];
                   }else if((del[i][j] > 0) && (del[i][j] > (box[j][j]/2))){
                       del[i][j] = del[i][j] - box[j][j];
                   }
                }
                
                msdX = msdX + (del[i][0]*del[i][0]);          msdY = msdY + (del[i][1]*del[i][1]);           msdZ = msdZ + (del[i][2]*del[i][2]);
                msd = msd + (del[i][0]*del[i][0]) + (del[i][1]*del[i][1]) + (del[i][2]*del[i][2]);
            }
            msdX /= nOxygen;  msdY /= nOxygen;  msdZ /= nOxygen;
    }
    private void write_xvg_mean_sd(float fc, int r, float[][] arr1, float[][] arr2, float[][] box, String output){
            // The mean and standard deviation for this method was tested and validated with calculations in R.
            int nOxygen = r/3;
            float[][] del = new float[nOxygen][3];
            float stdvX, stdvY, stdvZ;
            stdvX = stdvY = stdvZ = 0; 
            int i, j;
            for(i = 0; i < nOxygen; i++){
                del[i][0] = arr2[i*3][0] - arr1[i*3][0];    del[i][1] = arr2[i*3][1] - arr1[i*3][1];     del[i][2] = arr2[i*3][2] - arr1[i*3][2];
                if(output.equals("crd")){
                   for(j = 0; j < 3; j++){ 
                      if((del[i][j] < 0) && (del[i][j] < (-1*box[j][j]/2))){
                          del[i][j] = del[i][j] + box[j][j];
                      }else if((del[i][j] > 0) && (del[i][j] > (box[j][j]/2))){
                          del[i][j] = del[i][j] - box[j][j];
                      }
                   }
                }
                sumX = sumX + (del[i][0]);          sumY = sumY + (del[i][1]);           sumZ = sumZ + (del[i][2]);
            }
            sumX /= nOxygen;  sumY /= nOxygen;  sumZ /= nOxygen;
            for(i = 0; i < nOxygen; i++){
                stdvX = stdvX + ((del[i][0] - sumX)*(del[i][0] - sumX));
                stdvY = stdvY + ((del[i][1] - sumY)*(del[i][1] - sumY));
                stdvZ = stdvZ + ((del[i][2] - sumZ)*(del[i][2] - sumZ));
            }           
            stdvX = (float)Math.sqrt(stdvX/(nOxygen-1));  stdvY = (float)Math.sqrt(stdvY/(nOxygen-1));  stdvZ = (float)Math.sqrt(stdvZ/(nOxygen-1));
            if(output.equals("crd")){
               doutC.println(fc+"\t"+sumX+"\t"+stdvX+"\t"+sumY+"\t"+stdvY+"\t"+sumZ+"\t"+stdvZ);
            }else if(output.equals("vel")){
               doutV.println(fc+"\t"+sumX+"\t"+stdvX+"\t"+sumY+"\t"+stdvY+"\t"+sumZ+"\t"+stdvZ);
            }else if(output.equals("frc")){
               doutF.println(fc+"\t"+sumX+"\t"+stdvX+"\t"+sumY+"\t"+stdvY+"\t"+sumZ+"\t"+stdvZ);
            } 
    }
   // private float mean(int nO){return 0.0;} 
    private void read_header(String line, String[] frameData){
         try{
            String[] tempArr;
            String temp;
            StringTokenizer token;
            token = new StringTokenizer(line);
            int c, r, i, j, tkns;
            if(token.countTokens() != 3){ 
               System.out.println("ERROR:A line that does not have 3 tokens can not be the begining of a frame.");
            }else{
                temp = token.nextToken(); temp = token.nextToken();
                tempArr = token.nextToken().split(":");
                frameData[0] = tempArr[0];

                line = dFile.readLine();                   token = new StringTokenizer(line);
                temp = token.nextToken();                  frameData[1] = token.nextToken();
                temp = token.nextToken();                  frameData[2] = token.nextToken();
                tempArr = token.nextToken().split("=");    frameData[3] = tempArr[1];
                temp = token.nextToken();                  frameData[4] = token.nextToken();
            }
         }catch( IOException e ){System.err.println( e );}
    }
    private void read_frame(float[][] box, float[][] crd, float[][] vel, float[][] frc){
         try{
            String[] tempArr;
            String line, temp;
            StringTokenizer token;
            int c, r, i, j, tkns; 
            line = dFile.readLine();                    token = new StringTokenizer(line);
            if(token.nextToken().equals("box")){
               tempArr = token.nextToken().split("\\(|x|\\)|:");
               r = Integer.parseInt(tempArr[1]);    c = Integer.parseInt(tempArr[2]);
               read_Array(box, r);
            }
            line = dFile.readLine();                    token = new StringTokenizer(line);
            if(token.nextToken().equals("x")){
               tempArr = token.nextToken().split("\\(|x|\\)|:");
               r = Integer.parseInt(tempArr[1]);    c = Integer.parseInt(tempArr[2]);
               read_Array(crd, r);
            }
            line = dFile.readLine();                    token = new StringTokenizer(line);
            if(token.nextToken().equals("v")){
               tempArr = token.nextToken().split("\\(|x|\\)|:");
               r = Integer.parseInt(tempArr[1]);    c = Integer.parseInt(tempArr[2]);
               read_Array(vel, r);
            }
            line = dFile.readLine();                    token = new StringTokenizer(line);
            if(token.nextToken().equals("f")){
               tempArr = token.nextToken().split("\\(|x|\\)|:");
               r = Integer.parseInt(tempArr[1]);    c = Integer.parseInt(tempArr[2]);
               read_Array(frc, r);
            }
         }catch( IOException e ){System.err.println( e );}
    }
    private void read_Array(float[][] arr, int r){
         try{
            int i, j, tkns;
            String[] tempArr;
            String line, temp; 
            StringTokenizer  token;
            for(i = 0; i < r; i++){   //Turn this for loop into a function because it will be reused 3 times.
                line = dFile.readLine();                    token = new StringTokenizer(line);             tkns = token.countTokens();
                if(tkns == 4){
                   temp = token.nextToken();
                   tempArr = token.nextToken().split("\\]|=|\\{|,");
                   arr[i][0] = Float.parseFloat(tempArr[3]);
                   tempArr = token.nextToken().split(",");
                   arr[i][1] = Float.parseFloat(tempArr[0]);
                   tempArr = token.nextToken().split("\\}");
                   arr[i][2] = Float.parseFloat(tempArr[0]);
                }else if(tkns == 5){
                   temp = token.nextToken();        temp = token.nextToken();
                   tempArr = token.nextToken().split("\\]|=|\\{|,");
                   arr[i][0] = Float.parseFloat(tempArr[0]);
                   tempArr = token.nextToken().split(",");
                   arr[i][1] = Float.parseFloat(tempArr[0]);
                   tempArr = token.nextToken().split("\\}");
                   arr[i][2] = Float.parseFloat(tempArr[0]);
                }
            }
         }catch( IOException e ){System.err.println( e );}
    }
}
