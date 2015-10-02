package Aguan.files;
import Utilities.binaryManipulation;
import java.io.*;

public class trr extends trajectory {
    private boolean is_xtc_read, bigEndian;
    private int magic, step, atoms;
    private float time;
    private float[][] box;
    private float[][] coords, velocities, forces;
    private int[] frame_index;
    private int frame_counter, frame_atoms; 
    private int box_size, wrapper_size, header_size;
    private long estimated_file_length;
    private int bytes_in_one_frame, number_of_frames;
    /**
    * @param xtc The path to the xtc file.
    * Class constructed from file description in wed address:
    * http://manual.gromacs.org/online/xtc.html
    */
    public trr(String trr){
        super(trr);
        box = new float[3][3];
        frame_counter = 0;
        box_size = 36;
        wrapper_size = 4;
        header_size = 16;
    }
    public void trajquery(){
        read_frame();
        print_header();
        est_file_length();
        if((super.bytes_read_file_length % bytes_in_one_frame) == 0){
            System.out.println("Number of frames: "+number_of_frames);
            System.out.println("Estimated file length equal read file length: "+estimated_file_length+" == "+super.bytes_read_file_length);
        }else{
            System.out.println("Frame size and total file size are not multiples.");
            System.out.println(estimated_file_length+" != "+super.bytes_read_file_length);
        }
        super.closeInputStream();
    }
    public void mod_header(String xtc2Path, String field, int val){
        read_frame();
        System.out.println("Old header:");
        print_header();
              if(field.equals("magic")){   magic = val;
        }else if(field.equals("atoms")){   atoms = val;
        }else{ System.out.println("Only 'magic' and 'atoms' fields can be modified. 'Step' and 'time' vary for each frame.");}
        est_file_length();
        System.out.println("Modified DCD:");
        print_header();
        if((super.bytes_read_file_length % bytes_in_one_frame) == 0){
            gen_TRR_body_pointers();
            openOutputStream(xtc2Path);
            write_frame();
            for(int i = 0; i < number_of_frames; i++){
                read_frame();
                write_frame();
            }
            System.out.println("Number of frames: "+number_of_frames);
            System.out.println("Estimated file length equal read file length: "+estimated_file_length+" == "+super.bytes_read_file_length);
        }else{
            System.out.println("Frame size and total file size are not multiples.");
            System.out.println(estimated_file_length+" != "+super.bytes_read_file_length);
        }
        super.closeInputStream();
    }
    private void write_TRR(){
           try{
              write_frame();
              byte[] frame = new byte[(int)bytes_in_one_frame];
              int xRead;
              for(int y = 1; y < number_of_frames; y++){
                 xRead = dis.read(frame,0,frame.length);
                 outdcd.write(frame);
              }
           }catch( IOException e ){System.err.println( e );}
    }
    //outdcd.write(binaryManipulation.intToDWord_LE(NSET));
    private void gen_TRR_body_pointers(){
           frame_index = new int[number_of_frames];
           frame_index[0] = 0;
           System.out.println("0 "+frame_index[0]);
           for(int i = 1; i < number_of_frames; i++){
               frame_index[i] = frame_index[i-1] + bytes_in_one_frame;
               System.out.println(i+" "+frame_index[i]);
           }
    }
    private void read_frame(){
        try{
           int i, j;
  //         magic = dis.readInt();
  //         atoms = dis.readInt();
  //         step = dis.readInt();
  //         time = dis.readFloat(); HERE
           System.out.println("magic: "+dis.readInt()+", atoms: "+dis.readInt()+", step: "+dis.readInt()+", time: "+dis.readFloat());
           for(i = 0; i < 3; i++){
               System.out.print("|");
               for(j = 0; j < 3; j++){box[i][j] = dis.readFloat(); System.out.print(box[i][j]+", ");}
               System.out.println("|");
           }
           frame_atoms = dis.readInt();
           System.out.println("frmae_atoms: "+frame_atoms);
           coords = new float[frame_atoms][3];
           for(i = 0; i < frame_atoms; i++){
               System.out.print("|");
               for(j = 0; j < 3; j++){coords[i][j] = dis.readFloat(); System.out.print(coords[i][j]+", ");}
               System.out.println("|");
           }
           frame_counter++;
        }catch( IOException e ){System.err.println( e );}
    }
    private void write_frame(){
        try{
           int i, j;
           dout.writeInt(magic);
           dout.writeInt(atoms);
           dout.writeInt(step);
           dout.writeFloat(time);
           for(i = 0; i < 3; i++){
               for(j = 0; j < 3; j++){dout.writeFloat(box[i][j]);}
           }
           dout.writeInt(frame_atoms);
           for(i = 0; i < frame_atoms; i++){
               for(j = 0; j < 3; j++){dout.writeFloat(coords[i][j]);}
           }
        }catch( IOException e ){System.err.println( e );}
    }
    private void print_header(){
        System.out.println("TRAJQUERY:");
        System.out.println("magic: "+magic);
        System.out.println("atoms: "+atoms);
        System.out.println("step: "+step);
        System.out.println("time: "+time);
    }
    private void est_file_length(){
        bytes_in_one_frame = (((atoms*12)+box_size+wrapper_size+header_size));
        number_of_frames = (int)super.bytes_read_file_length/bytes_in_one_frame;
        estimated_file_length = bytes_in_one_frame*(number_of_frames);
    }
}
