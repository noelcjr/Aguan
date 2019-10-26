/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package Aguan.files;

import java.io.*;
import java.util.*;
import Aguan.files.crd;
/**
 *
 * @author Noel Carrascal
 */
public class checkCHARMM {
    private crd coord;
    private boolean allTestsPassed;
    public checkCHARMM(String[] args){
           System.out.println("Number of arguments "+args.length);
           if(args.length == 2){
              System.out.println("Red file "+args[1]);
              coord = new crd();
              coord.readCRD(args[1]); 
              allTestsPassed = false;
           }else{
              System.out.println("Number of parameters entered:"+args.length);
              System.out.println("One parameter needed after '-v'.");
              System.out.println("This check assumes that only water  molecules are");
              System.out.println("present in the CRD file.");
              System.out.println("   1- Path to CRD file.");
           }
    }
}
