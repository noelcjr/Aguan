/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package Aguan.files;

import java.io.*;
import java.util.*;
/**
 *
 * @author Noel Carrascal
 */
public class compareFiles {
    private FileReader FR;
    private BufferedReader dataFile;
    private String line;
    private boolean allTestsPassed;
    
    public compareFiles(String[] args){
        try{
           if(args.length == 2){
              allTestsPassed = true;
              FR =  new FileReader(args[1]);
              dataFile = new BufferedReader( FR );
              while( ( line = dataFile.readLine( ) ) != null ){
                  StringTokenizer token = new StringTokenizer(line);
                  cmpFls(token.nextToken(),token.nextToken());
              }
              System.out.println("+++++++++++++++++++++++++++++++++++");
              System.out.println("allTestsPassed?:"+allTestsPassed);
           }else{
              System.out.println("Following the -t copmmand, enter a path to a Test directory");
              System.out.println("where the testFiles.txt file is found.");
           }
        }catch( IOException e ){System.err.println( e );}
    }
    public void cmpFls(String file1, String file2){
        try{
           FileReader test, ref;
           BufferedReader d1, d2;
           test =  new FileReader(file1);
           d1 = new BufferedReader( test );
           ref =  new FileReader(file2);
           d2 = new BufferedReader( ref );
           String line1, line2;
           boolean passedTest = true;
           boolean readlines = true;
           //for(int i = 0; i < lines; i++){
           while(readlines){
               line1 = d1.readLine();
               line2 = d2.readLine();
               if((line1 != null) && (line2 != null)){
                   if(line1.equals(line2)){
                      
                   }else{
                      passedTest = false;
                      allTestsPassed = false;
//                      System.out.println("   Line "+line1);
//                      System.err.println("    and "+line2+" are unequal.");
                   }
               }else{
                   readlines = false;
               }
           }
           System.out.println("PassTest:"+passedTest+" vimdiff "+file1+" "+file2);
        }catch( IOException e ){System.err.println( e );}
    }
}
