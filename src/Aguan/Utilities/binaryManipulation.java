package Aguan.Utilities;
import java.io.*;

public class binaryManipulation{
    private binaryManipulation(){
        throw new AssertionError();
    }
    public static byte[] intToDWord_BE(int i) {
		byte[] dword = new byte[4];
		dword[0] = (byte) (i & 0x00FF);
		dword[1] = (byte) ((i >> 8) & 0x000000FF);
		dword[2] = (byte) ((i >> 16) & 0x000000FF);
		dword[3] = (byte) ((i >> 24) & 0x000000FF);
		return dword;
    }
    public static byte[] intToDWord_LE(int i) {
		byte[] dword = new byte[4];
		dword[0] = (byte) ((i >> 24) & 0x000000FF);
		dword[1] = (byte) ((i >> 16) & 0x000000FF);
		dword[2] = (byte) ((i >> 8) & 0x000000FF);
		dword[3] = (byte) (i & 0x00FF);
		return dword;
    }
    private static final int MASK = 0xff;  
    public static byte[] intToByteArray(int param) {
          byte[] result = new byte[4];
          for (int i = 0; i < 4; i++) {
              int offset = (result.length - 1 - i) * 4;
              result[i] = (byte) ((param >>> offset) & MASK);
          }
          return result;
    }
    public static int byteArrayToInt(byte test[]){
          int bits = 0;
          int MASK = 0xff;
          int i = 0;
          for(int shifter = 3; shifter >= 0; shifter--){
              bits |= ((int) test[i] & MASK) << (shifter * 8);
	      i++;
	  }
          return bits;
    }
    public static String byteArrayToString(byte[] byteArray) {
          StringBuilder sb = new StringBuilder("[");
          if(byteArray == null) {
              throw new IllegalArgumentException("byteArray must not be null");
          }
          int arrayLen = byteArray.length;
          for(int i = 0; i < arrayLen; i++) {
              sb.append(byteArray[i]);
              if(i == arrayLen - 1) {
                  sb.append("]");
              } else{
                  sb.append(", ");
              }
          }
          return sb.toString();
    }
}
