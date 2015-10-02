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
              System.out.println("Four parameters needed:");
              System.out.println("   1- Path to edr_out file.");
              System.out.println("   2- Path to trr_out file.");
              System.out.println("   3- Number of frames.");
           }
    }
}
