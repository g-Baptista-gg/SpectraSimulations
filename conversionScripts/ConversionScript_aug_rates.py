from datetime import datetime
import sys

Shelli = []
LowerConfigi = []
JJi = []
Eigeni = []
Configi = []
Percentagei = []

Shellf = []
LowerConfigf = []
JJf = []
Eigenf = []
Configf = []
Percentagef = []

#ALL VALUES IN eV WHEN APPLICABLE
Energies = []
Rate = []
TotRateIS = []
BranchingRatio = []

def main():
    if len(sys.argv) != 2:
        print("Expected charge state in input... terminating.")
        return
    
    chargeState = sys.argv[1]
    
    with open("rates_raw/" + chargeState + "/Fe_" + chargeState + "_rates_auger.txt", "r") as lines:
        header = lines.readline().strip() #header line
        #print(header.split("\t"))
        for i, line in enumerate(lines):
            values = line.strip().split("\t")
            #print(i)
            
            Shelli.append(values[1].strip())
            LowerConfigi.append(values[2].strip())
            JJi.append(int(values[3].strip()))
            Eigeni.append(values[4].strip())
            Configi.append(values[5].strip())
            try:
                Percentagei.append(float(values[6].strip()))
            except:
                Percentagei.append(float(0.0))
            
            if(len(values) == 17):
                Shellf.append(values[7].strip())
                LowerConfigf.append(values[8].strip())
                JJf.append(int(values[9].strip()))
                Eigenf.append(values[10].strip())
                Configf.append(values[11].strip())
                try:
                    Percentagef.append(float(values[12].strip()))
                except:
                    Percentagef.append(float(0.0))
                
                Energies.append(values[13].strip())
                Rate.append(values[14].strip())
                TotRateIS.append(values[15].strip())
                BranchingRatio.append(values[16].strip())
            else:#(len(values) == 15)
                Shellf.append(values[7].strip())
                LowerConfigf.append(values[8].strip())
                JJf.append(int(values[9].strip()))
                Eigenf.append(values[10].strip())
                Configf.append("")
                Percentagef.append(float(0.0))
                
                Energies.append(values[11].strip())
                Rate.append(values[12].strip())
                TotRateIS.append(values[13].strip())
                BranchingRatio.append(values[14].strip())
            
    for record in range(len(Shelli)):
        if Shelli[record] == "1s":
            Shelli[record] = "K1"
        
        if Shelli[record] == "2s":
            Shelli[record] = "L1"
        if "2s" in Shellf[record]:
            Shellf[record] = Shellf[record].replace("2s", "L1").replace("_", "")
        
        if Shelli[record] == "3s":
            Shelli[record] = "M1"
        if "3s" in Shellf[record]:
            Shellf[record] = Shellf[record].replace("3s", "M1").replace("_", "")
        
        if Shelli[record] == "4s":
            Shelli[record] = "N1"
        if "4s" in Shellf[record]:
            Shellf[record] = Shellf[record].replace("4s", "N1").replace("_", "")
        
        if Shelli[record] == "2p":
            if "2p*1" in Configi[record]:
                Shelli[record] = "L2"
            elif "2p*2" in Configi[record]:
                Shelli[record] = "L3"
        if "2p" in Shellf[record]:
            if "2p*1" in Configf[record] and "2p4" in Configf[record]:
                Shellf[record] = Shellf[record].replace("2p", "L2").replace("_", "")
            elif "2p*2" in Configf[record] and "2p3" in Configf[record]:
                Shellf[record] = Shellf[record].replace("2p", "L3").replace("_", "")
            elif "2p4" in Configf[record] and "2p*1" not in Configf[record]:
                Shellf[record] = Shellf[record].replace("2p", "L2").replace("_", "")
            elif "2p2" in Configf[record]:
                Shellf[record] = Shellf[record].replace("2p", "L3").replace("_", "")
            elif "2p*1" in Configf[record] and "2p3" in Configf[record]:
                Shellf[record] = "L2L3"
        
        if Shelli[record] == "3p":
            if "3p*1" in Configi[record]:
                Shelli[record] = "M2"
            elif "3p*2" in Configi[record]:
                Shelli[record] = "M3"
        if "3p" in Shellf[record]:
            if "3p*1" in Configf[record] and "3p4" in Configf[record]:
                Shellf[record] = Shellf[record].replace("3p", "M2").replace("_", "")
            elif "3p*2" in Configf[record] and "3p3" in Configf[record]:
                Shellf[record] = Shellf[record].replace("3p", "M3").replace("_", "")
            elif "3p4" in Configf[record] and "3p*1" not in Configf[record]:
                Shellf[record] = Shellf[record].replace("3p", "M2").replace("_", "")
            elif "3p2" in Configf[record]:
                Shellf[record] = Shellf[record].replace("3p", "M3").replace("_", "")
            elif "3p*1" in Configf[record] and "3p3" in Configf[record]:
                Shellf[record] = "M2M3"
        
        if Shelli[record] == "3d":
            if "3d*2" in Configi[record]:
                Shelli[record] = "M4"
            elif "3d*1" in Configi[record]:
                Shelli[record] = "M5"
        if "3d" in Shellf[record]:
            if "3d*2" in Configf[record] and "3d3" in Configf[record]:
                Shellf[record] = Shellf[record].replace("3d", "M4").replace("_", "")
            elif "3d*1" in Configf[record] and "3d4" in Configf[record]:
                Shellf[record] = Shellf[record].replace("3d", "M5").replace("_", "")
            elif "3d3" in Configf[record] and "3d*2" in Configf[record]:
                Shellf[record] = Shellf[record].replace("3d", "M4").replace("_", "")
            elif "3d4" in Configf[record]:
                Shellf[record] = Shellf[record].replace("3d", "M5").replace("_", "")
            elif "3d*3" in Configf[record] and "3d5" in Configf[record]:
                Shellf[record] = "M4M5"
        

    with open("rates_converted/" + chargeState + "/26-augrate_" + chargeState + ".out", "w") as output:
        output.write("# Atomic number Z= 26  Date:" + datetime.today().strftime('%d-%m-%Y') + "\n\n")
        output.write("# Register Shell IS\t   Configuration IS \t\t\t\t2JJ IS \t\t Eigi \t  Higher Config IS \t\t\t\t\t percentage IS ---> \t Shell FS \t  Configuration FS \t\t\t\t\t2JJ FS \t\t Eigf \t  Higher Config FS \t\t\t\t  percentage FS \t\tEnergy(eV)\t\t\t rate(s-1)\t\t\ttotal rate from IS\t\t\tBranchingRatio\n")
        
        for i in range(len(Shelli)):
            output.write(f'\t{str(i+1):<2}\t{Shelli[i]:>5}\t\t   {LowerConfigi[i]:<29}\t{str(JJi[i]):>5}\t{Eigeni[i]:>10}\t{Configi[i]:>35}\t{Percentagei[i]:>10} \t\t   ---> \t{Shellf[i]:>5}\t{LowerConfigf[i]:>35}\t{str(JJf[i]):>5}\t{Eigenf[i]:>10}\t{Configf[i]:>35}\t{Percentagef[i]:>10}\t\t\t\t{Energies[i]:<18}\t {Rate[i]:<19}\t{TotRateIS[i]:<25}\t{BranchingRatio[i]:<23}\n')

if __name__ == "__main__":
   main()