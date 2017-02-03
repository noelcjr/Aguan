/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package Aguan.files;

import java.io.*;
import java.util.*;
import Aguan.files.trr_out;
/**
 *
 * @author Noel Carrascal
 */
public class checkGromacs {
    private trr_out to;
    private boolean allTestsPassed;
    
    public checkGromacs(String[] args){
           if(args.length == 8){
              to = new trr_out(args);
              allTestsPassed = false;
           }else{
              System.out.println("Number of parameters entered:"+args.length);
              System.out.println("Seven parameters needed after '-v':");
              System.out.println("   1- Path to edr_out file.");
              System.out.println("   2- Path to trr_out file.");
              System.out.println("   3- Number of frames.");
              System.out.println("   4- Number of Atoms.");
              System.out.println("   5- Target Temperature.");
              System.out.println("   6- rCut for cutoff.");
              System.out.println("   7- Tau for Berendsen thermostat.");
           }
    }
}
