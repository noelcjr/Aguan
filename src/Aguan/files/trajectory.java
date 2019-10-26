package Aguan.files;

public class trajectory extends file {
    public String traj_type;
    public trajectory(String traj_file){
            super(traj_file);
            super.file_type = "trajectory";
    }
}
