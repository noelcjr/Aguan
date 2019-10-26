package Aguan.files;
import java.io.*;

public class file {
    private File fl, outF;
    public FileOutputStream outdcd;
    public DataOutputStream dout;
    private InputStream is;
    public DataInputStream dis;
    public PrintWriter PW;
    private FileWriter FW;
    private FileReader FR;
    private BufferedReader BR;
    private boolean is_file_read = false;
    private boolean os_file_write = false;
    public long bytes_read_file_length;
    public String file_type = "";
    public String path;
    //////////////////////////////////////////////
    public file(){}
    public file(String fl_path){
           path = fl_path;
           fl =  new File(fl_path);
           bytes_read_file_length = fl.length();
           is_file_read = true;
    }
    public void openInputStream(){
        try{
           is = new FileInputStream(path);
           dis = new DataInputStream(is);
        }catch( IOException e ){System.err.println( e );}
    }
    public void closeInputStream(){ 
        try{ 
           is.close(); 
           dis.close(); 
        }catch( IOException e ){System.err.println( e );}
    }
    public void openOutputStream(String file){
        try{
           outF = new File(file);
           outdcd = new FileOutputStream (file);
           dout = new DataOutputStream (outdcd);
           os_file_write = true;
        }catch( IOException e ){System.err.println( e );}
    }
    public void openInputWriter(String file_name){
//        try{
           
  //      }catch( IOException e ){System.err.println( e );}
    }
    public void openOutputWriter(String file_name){
        try{
           FW = new FileWriter(file_name);
           PW = new PrintWriter(FW,true);
        }catch( IOException e ){System.err.println( e );}
    }
}
