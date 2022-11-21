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
Width = []
MultipoleNum = []
TotRateIS = []
BranchingRatio = []
MultipoleRates = []

h = 4.135667696 * 10**(-15)

def main():
    with open("Cu_rates_satellites.txt", "r") as lines:
        header = lines.readline().strip() #header line
        #print(header.split("\t"))
        for i, line in enumerate(lines):
            values = line.strip().split("\t")
            
            Shelli.append(values[1].strip())
            LowerConfigi.append(values[2].strip())
            JJi.append(int(values[3].strip()))
            Eigeni.append(values[4].strip())
            Configi.append(values[5].strip())
            try:
                Percentagei.append(float(values[6].strip()))
            except:
                Percentagei.append(float(0.0))
            
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
            Width.append(float(Rate[-1]) * h)
            MultipoleNum.append(values[15].strip())
            TotRateIS.append(values[16].strip())
            BranchingRatio.append(values[17].strip())
            MultipoleRates.append("\t".join(values[18:]).strip())
            
    for record in range(len(Shelli)):
        if "1s" in Shelli[record]:
            Shelli[record] = Shelli[record].replace("1s", "K1").replace("_", "")
        if "1s" in Shellf[record]:
            Shellf[record] = Shellf[record].replace("1s", "K1").replace("_", "")
        
        if "2s" in Shelli[record]:
            Shelli[record] = Shelli[record].replace("2s", "L1").replace("_", "")
        if "2s" in Shellf[record]:
            Shellf[record] = Shellf[record].replace("2s", "L1").replace("_", "")
        
        if "3s" in Shelli[record]:
            Shelli[record] = Shelli[record].replace("3s", "M1").replace("_", "")
        if "3s" in Shellf[record]:
            Shellf[record] = Shellf[record].replace("3s", "M1").replace("_", "")
        
        if "4s" in Shelli[record]:
            Shelli[record] = Shelli[record].replace("4s", "N1").replace("_", "")
        if "4s" in Shellf[record]:
            Shellf[record] = Shellf[record].replace("4s", "N1").replace("_", "")
        
        if "2p" in Shelli[record]:
            if "2p*1" in Configi[record] and "2p4" in Configi[record]:
                Shelli[record] = Shelli[record].replace("2p", "L2").replace("_", "")
            elif "2p*2" in Configi[record] and "2p3" in Configi[record]:
                Shelli[record] = Shelli[record].replace("2p", "L3").replace("_", "")
            elif "2p4" in Configi[record] and "2p*1" not in Configi[record]:
                Shelli[record] = Shelli[record].replace("2p", "L2").replace("_", "")
            elif "2p2" in Configi[record]:
                Shelli[record] = Shelli[record].replace("2p", "L3").replace("_", "")
            elif "2p*1" in Configi[record] and "2p3" in Configi[record]:
                Shelli[record] = "L2L3"
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
        
        if "3p" in Shelli[record]:
            if "3p*1" in Configi[record] and "3p4" in Configi[record]:
                Shelli[record] = Shelli[record].replace("3p", "M2").replace("_", "")
            elif "3p*2" in Configi[record] and "3p3" in Configi[record]:
                Shelli[record] = Shelli[record].replace("3p", "M3").replace("_", "")
            elif "3p4" in Configi[record] and "3p*1" not in Configi[record]:
                Shelli[record] = Shelli[record].replace("3p", "M2").replace("_", "")
            elif "3p2" in Configi[record]:
                Shelli[record] = Shelli[record].replace("3p", "M3").replace("_", "")
            elif "3p*1" in Configi[record] and "3p3" in Configi[record]:
                Shelli[record] = "M2M3"
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
        
        if "3d" in Shelli[record]:
            if "3d*3" in Configi[record] and "3d6" in Configi[record]:
                Shelli[record] = Shelli[record].replace("3d", "M4").replace("_", "")
            elif "3d*4" in Configi[record] and "3d5" in Configi[record]:
                Shelli[record] = Shelli[record].replace("3d", "M5").replace("_", "")
            elif "3d6" in Configi[record] and "3d*2" in Configi[record]:
                Shelli[record] = Shelli[record].replace("3d", "M4").replace("_", "")
            elif "3d4" in Configi[record]:
                Shelli[record] = Shelli[record].replace("3d", "M5").replace("_", "")
            elif "3d*3" in Configi[record] and "3d5" in Configi[record]:
                Shelli[record] = "M4M5"
        if "3d" in Shellf[record]:
            if "3d*3" in Configf[record] and "3d6" in Configf[record]:
                Shellf[record] = Shellf[record].replace("3d", "M4").replace("_", "")
            elif "3d*4" in Configf[record] and "3d5" in Configf[record]:
                Shellf[record] = Shellf[record].replace("3d", "M5").replace("_", "")
            elif "3d6" in Configf[record] and "3d*2" in Configf[record]:
                Shellf[record] = Shellf[record].replace("3d", "M4").replace("_", "")
            elif "3d4" in Configf[record]:
                Shellf[record] = Shellf[record].replace("3d", "M5").replace("_", "")
            elif "3d*3" in Configf[record] and "3d5" in Configf[record]:
                Shellf[record] = "M4M5"

    with open("29-satrate.out", "w") as output:
        output.write("# Atomic number Z= 29  Date:" + datetime.today().strftime('%d-%m-%Y') + "\n\n")
        output.write("# Register Shelli\t2Ji\t\tEigi\t ---> \t Shellf\t\t2Jf\t\tEigf\tEnergy(eV)\t\t\tRate(s-1)\t\tWidth(eV)\n")
        
        for i in range(len(Shelli)):
            output.write(f'\t{str(i+1):<3}\t{Shelli[i]:>5}\t{str(JJi[i]):>5}\t{Eigeni[i]:>6}\t\t ---> \t{Shellf[i]:>5}\t{str(JJf[i]):>5}\t{Eigenf[i]:>6}\t\t{Energies[i]:<18}\t{Rate[i]:<15}\t{Width[i]:<15}\n')

if __name__ == "__main__":
   main()