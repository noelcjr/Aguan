package Aguan.files;
import Aguan.Utilities.binaryManipulation;
import java.io.*;

public class dcd extends trajectory {
    private String header;
    private int title_size, title_size2, first_int, second_int, Nindex, freeIndex, freeIndex2;
    private int N, NSET, ISTART, NSAVC, LENGTH, NAMNF, NTITLE, countBytes, frame_size;
    private float DELTA;
    private String[] ntitle;
    private int[] frame_index, freeIndexes, infoArr5, infoArr10;
    private boolean is_dcd_read, bigEndian;
    private long bytes_in_atom_coordinates, bytes_in_wrappers, bytes_in_full_header; 
    private long bytes_in_est_file_length;
    /**
    * @param dcd The path to the dcd file.
    */
    public dcd(String dcd){
        super(dcd);
        infoArr5 = new int[5]; infoArr10 = new int[10];
    }
    public void trajquery(){
        System.out.printf("DCD file %s \n",read_dcd());
        print_header();
        closeInputStream();
    }
    public void mod_header(String dcd2Path, String field, int val){
        System.out.printf("DCD file %s \n",read_dcd());
        print_header();
        mod_field(field, val);
        est_file_length();
        print_header();
        if(bytes_in_est_file_length == super.bytes_read_file_length){
           gen_DCD_body_pointers();
           openOutputStream(dcd2Path);
           write_header();
           write_body();
           System.out.println("Header modification fixed dcd parameters");
        }else{
           System.out.println("Header modification failed. File size and estimated file size");
           System.out.println("do not match.");
           System.out.println(bytes_in_est_file_length+" != "+super.bytes_read_file_length);
        }
        closeInputStream();
    }
    public void mod_header(String dcd2Path, String field, int val, String field2, int val2){
        System.out.printf("DCD file %s \n",read_dcd());
        print_header();
        mod_field(field, val);
        mod_field(field2, val2);
        est_file_length();
        print_header();
        if(bytes_in_est_file_length == super.bytes_read_file_length){
           gen_DCD_body_pointers();
           openOutputStream(dcd2Path);
           write_header();
           write_body();
           System.out.println("Header modification fixed dcd parameters");
        }else{
           System.out.println("Header modification failed. File size and estimated file size");
           System.out.println("do not match.");
           System.out.println(bytes_in_est_file_length+" != "+super.bytes_read_file_length);
        }
        closeInputStream();
    }
    public void write_header(){
           try{
              super.outdcd.write(binaryManipulation.intToDWord_LE(first_int));
              byte[] w = {67,79,82,68};
              outdcd.write(w);
              outdcd.write(binaryManipulation.intToDWord_LE(NSET));
              outdcd.write(binaryManipulation.intToDWord_LE(ISTART));
              outdcd.write(binaryManipulation.intToDWord_LE(NSAVC));
              int i;
              for(i = 0; i < 5; i++){outdcd.write(binaryManipulation.intToDWord_LE(infoArr5[i]));}
              outdcd.write(binaryManipulation.intToDWord_LE(NAMNF));
              outdcd.write(binaryManipulation.intToDWord_LE(1026003170));
              for(i = 0; i < 10; i++){outdcd.write(binaryManipulation.intToDWord_LE(infoArr10[i]));}
              outdcd.write(binaryManipulation.intToDWord_LE(second_int));
              outdcd.write(binaryManipulation.intToDWord_LE(title_size));
              outdcd.write(binaryManipulation.intToDWord_LE(NTITLE));
              for(i = 0; i < NTITLE; i++){outdcd.write(ntitle[i].getBytes());}
              outdcd.write(binaryManipulation.intToDWord_LE(title_size2));
              outdcd.write(binaryManipulation.intToDWord_LE(4));
              outdcd.write(binaryManipulation.intToDWord_LE(N));
              outdcd.write(binaryManipulation.intToDWord_LE(4));
              if (NAMNF != 0){
                 outdcd.write(binaryManipulation.intToDWord_LE(freeIndex));
                 for(i = 0; i< ((N)-(NAMNF)); ++i){outdcd.write(binaryManipulation.intToDWord_LE(freeIndexes[i]));}
                 outdcd.write(binaryManipulation.intToDWord_LE(freeIndex2));
              }
           }catch( IOException e ){System.err.println( e );}
    }
    public int read_int(){
           int temp = 0;
           try{
              temp = dis.readInt();
           }catch( IOException e ){System.err.println( e );}
           return temp;
    }
    public float read_float(){
           float temp = 0;
           try{
              temp = dis.readFloat();
           }catch( IOException e ){System.err.println( e );}
           return temp;
    }
    public byte[] get_frame(){
           byte[] frameXYZ = new byte[(int)frame_size];
           try{
              int xRead = dis.read(frameXYZ,0,frameXYZ.length);
           }catch( IOException e ){System.err.println( e );}
           return frameXYZ;
    }
    public void write_body(){
           try{
              byte[] frameX = new byte[(int)frame_size];
              int xRead;
              for(int y = 0; y < NSET; y++){
                 xRead = dis.read(frameX,0,frameX.length);
                 outdcd.write(frameX);
              }
           }catch( IOException e ){System.err.println( e );}
    }
    public void gen_DCD_body_pointers(){
           frame_index = new int[NSET];
           frame_index[0] = (int)bytes_in_full_header;
           frame_size = (N*12) + 24;
           int byteCount = frame_index[0];
           for(int i = 1; i < NSET; i++){
               frame_index[i] = frame_index[i-1] + frame_size;
           }
    }
    public String read_dcd(){
        try{
           int four, four2;
           first_int = dis.readInt();  
           System.out.println(countBytes+" "+first_int);    countBytes += 4;
           if(first_int != 84){
	      System.out.println("The first byte is not equal to 84. It is equal to "+first_int);
	      System.out.println("This means that the DCD file is eithe little endian or corrupted");
	      return "BAD DCD FORMAT: First Int is not 84.";
           }
           StringBuffer buffer = new StringBuffer();
           buffer.append((char)dis.readByte());   // 67
           buffer.append((char)dis.readByte());   // 79
           buffer.append((char)dis.readByte());   // 82
           buffer.append((char)dis.readByte());   // 68
       //    System.out.println("Int = "+dis.readInt());   // 1129271876
           header = buffer.toString();                      countBytes += 4;
           if(!header.equals("CORD")){return "BAD DCD FORMAT: Not CORD header.";}
           NSET = dis.readInt();                            countBytes += 4;
           ISTART = dis.readInt();                          countBytes += 4;
           NSAVC = dis.readInt();                           countBytes += 4;
           int i = 0;    int j = 0;
           for(i = 0; i < 5; i++){ infoArr5[i] = dis.readInt(); countBytes += 4;}
           NAMNF = dis.readInt();                           countBytes += 4;
           DELTA = dis.readFloat();                         countBytes += 4;
       //    System.out.println("DELTA + "+dis.readInt());   // 1026003170
           for(i = 0; i < 10; i++){
	       infoArr10[i] = dis.readInt();                countBytes += 4;
           }
           second_int = dis.readInt();                      countBytes += 4;
           if(second_int != 84){return "BAD DCD FORMAT: Second Int is not 84";}
           title_size = dis.readInt();                      countBytes += 4;
           if(((title_size-4)%80) == 0){
              NTITLE = dis.readInt();                       countBytes += 4;
              ntitle = new String[NTITLE];
              for(i = 0; i < NTITLE; i++){
                 StringBuffer buffer2 = new StringBuffer();
                 for(j = 0; j < 80; j++){
                     buffer2.append((char)dis.readByte());
                     countBytes++;
                 }
                 ntitle[i] = buffer2.toString();
              }
           }else{return "BAD DCD FORMAT: ((Title Size - 4) % 80) != 0";}
           title_size2 = dis.readInt();                      countBytes += 4;
           if(title_size != title_size2){return "BAD DCD FORMAT: Title_size 1 and do not match.";}
           four = dis.readInt();                             countBytes += 4;
           if(four != 4){return "BAD DCD FORMAT: Int after title not equal to 4 1.";}
           Nindex = countBytes;
           N = dis.readInt();                                countBytes += 4;
           four2 = dis.readInt();                            countBytes += 4;
           if(four2 != 4){return "BAD DCD FORMAT: Int after title not equal to 4 2.";}           
           if (NAMNF != 0){
               freeIndexes = new int[(N)-(NAMNF)];
               freeIndex = dis.readInt();                    countBytes += 4;
               if (freeIndex != ((N)-(NAMNF))*4){return "BAD DCD FORMAT: Inside free indexes1.";}
               for (i=0; i<((N)-(NAMNF)); ++i){
                   freeIndexes[i] = dis.readInt();           countBytes += 4;
               }
               freeIndex2 = dis.readInt();                   countBytes += 4;
               if(freeIndex2 != ((N)-(NAMNF))*4){return "BAD DCD FORMAT: Inside free indexes2.";}
           }
        }catch( IOException e ){System.err.println( e );}
        est_file_length();
        if(super.bytes_read_file_length == bytes_in_est_file_length){
           is_dcd_read = true;
           return "DCD file has been read without errors.";
        }else{
           is_dcd_read = false;
           return "ERROR: In read_dcd(). DCD estimated file size != real file size.";
        }
    }
    public void print_header(){
        System.out.println("TRAJQUERY:");
        System.out.println("NSET = "+NSET);
        System.out.println("ISTART = "+ISTART);
        System.out.println("N = "+N);
        System.out.println("NSAVC = "+NSAVC);
        System.out.println("NAMNF = "+NAMNF);
        System.out.println("DELTA = "+DELTA);
        System.out.println("NTITLE = "+NTITLE);
        System.out.println("title_size = "+title_size);
        int i;
        System.out.print("infoArr05: ");for(i = 0; i < infoArr5.length; i++){System.out.print(infoArr5[i]+" ");}    System.out.println();
        System.out.print("infoArr10: ");for(i = 0; i < infoArr10.length; i++){System.out.print(infoArr10[i]+" ");}  System.out.println();
        System.out.print("ntitle   : ");for(i = 0; i < NTITLE; i++){System.out.print(i+" "+ntitle[i]);}             System.out.println();
        System.out.println("Number of frames NSET = "+NSET);
        System.out.println("Number of atoms N = "+N);
        System.out.println("Bytes in full header                      = "+bytes_in_full_header);	
        System.out.println("Bytes used for holding Atoms coordinates  = "+bytes_in_atom_coordinates);
        System.out.println("Bytes used for holding wrappers           = "+bytes_in_wrappers);
        // file length =  header_length + Atoms_length + wrappers_legth
        System.out.println("    Estimated file Length                 = "+bytes_in_est_file_length);
        if(super.bytes_read_file_length != bytes_in_est_file_length){
           System.out.println("ERROR: In read_dcd(). DCD estimated file size != real file size.");
        }else{
           System.out.println("DCD estimated file size matches real size.");
        }
    }
    public void mod_field(String field, int val){
              if(field.equals("N")){        N = val;
        }else if(field.equals("NSET")){     NSET = val;
        }else if(field.equals("ISTART")){   ISTART = val;
        }else if(field.equals("NSAVC")){    NSAVC = val;
        }else if(field.equals("LENGTH")){   LENGTH = val;
        }else if(field.equals("NAMNF")){    NAMNF = val;
        }else if(field.equals("NTITLE")){   NTITLE = val;
        }else{ System.out.println("Wrong field to be modified."); }
    }
    public void est_file_length(){
        bytes_in_atom_coordinates = ((long)N*3*4*(long)NSET);
        bytes_in_wrappers = (long)((long)NSET*3*2*4);
        bytes_in_full_header = (long)(112+title_size);
        bytes_in_est_file_length = bytes_in_atom_coordinates + bytes_in_wrappers + bytes_in_full_header;
    }
    public int get_NSET(){
        return NSET;
    }
    public int get_N(){
        return N;
    }
}
