package Aguan.files;

public class coordinate extends file {
    public String coordinate_type;
    public coordinate(String coordinate_file){
           super(coordinate_file);
           super.file_type = "coordinate";
    }
}
