package Aguan.files;
import Aguan.TheMatrix.TheMatrix;
import java.io.*;

public class output extends file {
    public PrintWriter pout;
    public output(String output_file){
         super(output_file);
         super.file_type = "output";
         super.openOutputWriter(output_file);
         pout = super.PW;
    }
    public void PrintSummary(TheMatrix TM, double timeNow){
         pout.printf("%.4f %5.4f ",(timeNow*1.6),(TM.totEnergyVal*TM.ep[0]));
         pout.printf("%5.4f %5.4f %5.4f %5.4f %5.4f %5.4f ",(TM.potEnergyVal*TM.ep[0]),(TM.potEnergyVdwVal*TM.ep[0]),
                                                            (TM.potEnergyEeVal*TM.ep[0]),(TM.potEnergyRf1Val*TM.ep[0]),
                                                            (TM.potEnergyRf2Val*TM.ep[0]),(TM.potEnergyEfVal*TM.ep[0]));
         pout.printf("%5.4f %5.4f %5.4f ",(TM.kinEnergyVal*TM.ep[0]),(TM.rotKinEnergyVal*TM.ep[0]),(TM.trzKinEnergyVal*TM.ep[0]));
         pout.printf("%5.4f %5.4f %5.4f ",(TM.instantTemperature*78.42),(TM.translationalTemperature*78.42),(TM.rotationalTemperature*78.42));
         pout.printf("%5.4f \n",TM.pressure*TM.pre_convFact[0]);//mega pascal convertion factor
    }
    public void PrintVelocitiesOne(TheMatrix TM, double timeNow){
         pout.printf("%.4f ",(timeNow*1.6)); // Conversion factor from reduce units of time to femtoseconds
         pout.printf("%.4f %.4f %.4f ",TM.rvx[0],TM.rvy[0],TM.rvz[0]);
         pout.printf("%.8f %.8f %.8f\n",TM.wvx[0],TM.wvy[0],TM.wvz[0]);
    }
    public void PrintSummary_ave(TheMatrix TM, double timeNow){
        pout.printf("%.4f %5.4f ",(timeNow*1.6),(TM.totEnergySum*TM.ep[0]));
        pout.printf("%5.4f ",(TM.potEnergySum*TM.ep[0]));
        pout.printf("%5.4f ",(TM.potEnergyVdwSum*TM.ep[0]));
        pout.printf("%5.4f ",(TM.potEnergyEeSum*TM.ep[0]));
        pout.printf("%5.4f ",(TM.potEnergyRf1Sum*TM.ep[0]));
        pout.printf("%5.4f ",(TM.potEnergyRf2Sum*TM.ep[0]));
        pout.printf("%5.4f ",(TM.potEnergyEfSum*TM.ep[0]));
        pout.printf("%5.4f ",(TM.kinEnergySum*TM.ep[0]));
        pout.printf("%5.4f ",(TM.rotKinEnergySum*TM.ep[0]));
        pout.printf("%5.4f ",(TM.trzKinEnergySum*TM.ep[0]));
        pout.printf("%5.4f ",TM.aveTemperature*78.42);
        pout.printf("%5.4f ",(TM.trzTempAve*78.42));
        pout.printf("%5.4f ",(TM.rotTempAve*78.42));
        pout.printf("%5.4f \n",TM.pressureSum*TM.pre_convFact[0]);
        pout.flush();
    }
}
