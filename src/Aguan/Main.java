/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package Aguan;
import Aguan.files.input;
import Aguan.TheMatrix.TheMatrix;
import Aguan.files.restart;
import Aguan.Methods.rigidLeapFrogBerenTV;
import Aguan.Methods.minimization;
import Aguan.ntraj.generateMolecularOrientationTest;
import Aguan.Methods.pr_08_4;
import Aguan.files.compareFiles;
import Aguan.files.checkCHARMM;
import Aguan.ntraj.ntraj;
import Aguan.parameters;
import java.util.*;
/**
 *
 * @author Noel Carrascal
 */
public class Main {
    private restart RR;
    private rigidLeapFrogBerenTV RLFBTV;           private minimization min;
    private compareFiles cf;                       private ntraj ntj;
    private checkCHARMM cC;                        private pr_08_4 p84;
    public Main(String[] args){
                 if(args[0].equals("-s")){
                    RLFBTV = new rigidLeapFrogBerenTV(args);
           }else if(args[0].equals("-m")){
                    min = new minimization(args);
           }else if(args[0].equals("-o")){
                    //tt = new generateMolecularOrientationTest(TM);
           }else if(args[0].equals("-c")){
                    RR = new restart(args); 
           }else if(args[0].equals("-t")){
                    cf = new compareFiles(args); 
           }else if(args[0].equals("-n")){
                    ntj = new ntraj(args);
           }else if(args[0].equals("-v")){
                    cC = new checkCHARMM(args);
           }else if(args[0].equals("-z")){
                    p84 = new pr_08_4(args);
           }else{
              System.out.println("Type any of the following options after -n, -g or -t for intructions on how to use them.");
           }
    }
    public static void main(String[] args){
           if(args.length > 0){
              if(args[0].equals("-h")){
                 System.out.println("The first argument can only be:");
                 System.out.println("    -s  : This is for doing simulations. it requires two more paramater.");
                 System.out.println("               1) An input file. For the format of input read the manual.");
                 System.out.println("               2) A path to a blank output file.");
                 System.out.println("    -m  : Minimization. Takes a rst file and returns an rst file.");
                 System.out.println("    -n  : For ntraj trajectory analysis. if typed by itself it gives");
                 System.out.println("          instructions on how to use the trajectory analysis program.");
                 System.out.println("    -c  : For creating rst files automatically.");
                 System.out.println("    -t  : for Devoleping testing purposes. Users should not use this argument");
                 System.out.println("          because it requires a special directory to exist with a special");
                 System.out.println("          input file, or else it will generate an exception handling error.");
                 System.out.println("    -v  : Validations of trajectories, velocities and forces of my program to Gromacs results.");
                 System.out.println("          From trajectories obtained from Gromacs simulations, energies are calculated and matched");
                 System.out.println("          to energies from Gromacs outputs.");
                 System.out.println("    -o  : Molecular oreintation test. Function implemented, but marked for merging and modification.");
                 System.out.println("    -z  : Perform the same functionality as pr_08_4.c in the Rapaport text book. Used for testing.");
                 System.out.println("          It outputs rapaport's output. It requires two paramenters:");
                 System.out.println("             1) An imput file identical to the one needed for pr_08_4.c. Do not used other than pr_08_4.in.");
                 System.out.println("             2) A restart file generated from pr_08_4 modify to output coordinates, velocities,");
                 System.out.println("                accelerations, rotation matrices, angular velocity and angular accelerations");
              }else{
                 Main m = new Main(args);
              }
           }else{
              System.out.println("Type -h after \'java -jar my.jar\' for further instructions.");
           }
    }
}
