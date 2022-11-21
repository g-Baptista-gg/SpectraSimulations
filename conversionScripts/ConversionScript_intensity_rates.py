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
    with open("Cu_rates_radiative.txt", "r") as lines:
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
        if Shelli[record] == "1s":
            Shelli[record] = "K1"
        if Shellf[record] == "1s":
            Shellf[record] = "K1"
        
        if Shelli[record] == "2s":
            Shelli[record] = "L1"
        if Shellf[record] == "2s":
            Shellf[record] = "L1"
        
        if Shelli[record] == "3s":
            Shelli[record] = "M1"
        if Shellf[record] == "3s":
            Shellf[record] = "M1"
        
        if Shelli[record] == "4s":
            Shelli[record] = "N1"
        if Shellf[record] == "4s":
            Shellf[record] = "N1"
        
        if Shelli[record] == "2p":
            if "2p4" in Configi[record]:
                Shelli[record] = "L2"
            elif "2p*2" in Configi[record]:
                Shelli[record] = "L3"
        if Shellf[record] == "2p":
            if "2p4" in Configf[record]:
                Shellf[record] = "L2"
            elif "2p*2" in Configf[record]:
                Shellf[record] = "L3"
        
        if Shelli[record] == "3p":
            if "3p4" in Configi[record]:
                Shelli[record] = "M2"
            elif "3p3" in Configi[record] or "3p2" in Configi[record]:
                Shelli[record] = "M3"
        if Shellf[record] == "3p":
            if "3p4" in Configf[record]:
                Shellf[record] = "M2"
            elif "3p3" in Configf[record] or "3p2" in Configf[record]:
                Shellf[record] = "M3"
        
        if Shelli[record] == "3d":
            if "3d*3" in Configi[record]:
                Shelli[record] = "M4"
            elif "3d*4" in Configi[record]:
                Shelli[record] = "M5"
        if Shellf[record] == "3d":
            if "3d*3" in Configf[record]:
                Shellf[record] = "M4"
            elif "3d*4" in Configf[record]:
                Shellf[record] = "M5"

    with open("29-radrate.out", "w") as output:
        output.write("# Atomic number Z= 29  Date:" + datetime.today().strftime('%d-%m-%Y') + "\n\n")
        output.write("# Register Shelli\t2Ji\t\tEigi\t ---> \t Shellf\t\t2Jf\t\tEigf\tEnergy(eV)\t\t\tRate(s-1)\t\tWidth(eV)\n")
        
        for i in range(len(Shelli)):
            output.write(f'\t{str(i+1):<3}\t{Shelli[i]:>5}\t{str(JJi[i]):>5}\t{Eigeni[i]:>6}\t\t ---> \t{Shellf[i]:>5}\t{str(JJf[i]):>5}\t{Eigenf[i]:>6}\t\t{Energies[i]:<18}\t{Rate[i]:<15}\t{Width[i]:<15}\n')

if __name__ == "__main__":
   main()