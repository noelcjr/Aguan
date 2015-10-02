package Aguan.ntraj;
import java.io.*;
import Aguan.files.*;
//import Aguan.ntraj.*;
public class ntraj {
    private dcd DCD;
    private xtc XTC;
    private trr TRR;
    private diffusion CVF; 
    private String trj_typ, trj_typ2;
    public ntraj(String[] args){
           if(args.length > 1){
              trj_typ = args[2].substring((args[2].length()-3),args[2].length());
              if(args[1].equals("trajquery")){
                 if(trj_typ.equals("dcd")){
                    DCD = new dcd(args[2]);
                    DCD.openInputStream();
                    trajqueryDCD();
                 }else if(trj_typ.equals("xtc")){
                    XTC = new xtc(args[2]);
                    trajqueryXTC();
                 }else if(trj_typ.equals("trr")){
                    TRR = new trr(args[2]);
                    trajqueryTRR();
                 }else{
                    System.out.println("ERROR: Wrong ntraj file type for trajquery.");
                 }
              }else if(args[1].equals("diffusion")){
                 if(trj_typ.equals("dcd")){
                    CVF = new diffusion(args[2],"dcd");
                    CVF.read_dcd_file();
                 }else if(trj_typ.equals("txt")){
                    trj_typ = args[2].substring((args[2].length()-7),args[2].length());
                    if(trj_typ.equals("trr.txt")){
                       CVF = new diffusion(args[2],"trr");
                       CVF.read_trr_file();
                    }else{
                       System.out.println("File does not end with trr.txt. Make sure a trr.txt is generated from gmxdump -f *trr > *trr.txt");
                    }
                 }else{
                    System.out.println("ERROR: Wrong ntraj file type for diffusion.");
                 }              
              }else if(args[1].equals("mod_header")){
                 trj_typ2 = args[3].substring((args[3].length()-3),args[3].length());
                 if(trj_typ.equals(trj_typ2)){
                    if(trj_typ.equals("dcd")){
                       DCD = new dcd(args[2]);
                       DCD.openInputStream(); 
                       if(args.length == 6){
                          mod_headerDCD(args[3], args[4], Integer.parseInt(args[5]));
                       }else if(args.length == 8){
                          System.out.println("inside ntraj before mod_header is called.");
                          mod_headerDCD(args[3], args[4], Integer.parseInt(args[5]), args[6], Integer.parseInt(args[7]));
                       }
                    }else if(trj_typ.equals("xtc")){
                       XTC = new xtc(args[2]);
                       mod_headerXTC(args[3], args[4], Integer.parseInt(args[5]));
                    }else if(trj_typ.equals("trr")){
                       TRR = new trr(args[2]);
                       mod_headerTRR(args[3], args[4], Integer.parseInt(args[5]));
                    }else{   System.out.println("ERROR: Wrong ntraj file type.");}
                 }else{
                    System.out.println("ERROR: Traj files' types are different");
                 }
              }else{
                   System.out.println("Wrong commands following -n. Type -n only for copmmand options.");
              }
           }else{
              System.out.println("ntraj commands that need to follow -n:");
              command_description();
           }
    }
    private void trajqueryDCD(){                               DCD.trajquery();}
    private void trajqueryXTC(){                               XTC.trajquery();}
    private void trajqueryTRR(){                               TRR.trajquery();}
    private void mod_headerDCD(String trj2, String o, int val){DCD.mod_header(trj2,o,val);}
    private void mod_headerDCD(String trj2, String o, int val, String o2, int val2){DCD.mod_header(trj2,o,val,o2,val2);}
    private void mod_headerXTC(String trj2, String o, int val){System.out.println("mod_header for XTC files not implemented.");
                                                               System.out.println("Compression algorithm and byte count not known yet.");}//XTC.mod_header(trj2,o,val);}
    private void mod_headerTRR(String trj2, String o, int val){System.out.println("mod_header for TRR files not implemented.");}//TRR.mod_header(trj2,o,val);}
    private void command_description(){
            System.out.println("    trajquery  and a path to a dcd, xtc or trr file.");
            System.out.println("    diffusion  and a path to a file that end with trr.txt and that is obtained from gmxdump -f *.trr command.");
            System.out.println("    mod_header followed by a path to a dcd (xtc and trr files not supported yet) file, a field name to be ");
            System.out.println("               modified(N,NSET,ISTART,NSAVC,LENGTH,NAMNF,NTITLE), and the value for the field.");
    }
}
