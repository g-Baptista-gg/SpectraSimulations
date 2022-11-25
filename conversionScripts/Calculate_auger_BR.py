import sys

RadShelli = []
RadJJi = []
RadShellf = []
RadJJf = []
RadRate = []

AugShelli = []
AugJJi = []
AugEigvi = []
AugShellf = []
AugJJf = []
AugEigvf = []
AugRate = []

Shells = []
ShellsA = []
RadTotRate = []
AugTotRate = []
BR = []

Configs = []

def main():
    with open("29-radrate.out", "r") as intFile:
        intFile.readline()
        intFile.readline()
        intFile.readline()
        
        for line in intFile:
            values = line.strip().split()
            values = list(filter(None, values))
            
            RadShelli.append(values[1].strip())
            RadJJi.append(values[2].strip())
            
            RadShellf.append(values[5].strip())
            RadJJf.append(values[6].strip())
            
            RadRate.append(values[9].strip())
    
    with open("29-augerrate.out", "r") as augFile:
        augFile.readline()
        augFile.readline()
        augFile.readline()
        
        for line in augFile:
            values = line.strip().split()
            values = list(filter(None, values))
            
            AugShelli.append(values[1].strip())
            AugJJi.append(values[2].strip())
            AugEigvi.append(values[3].strip())
            
            AugShellf.append(values[5].strip())
            AugJJf.append(values[6].strip())
            AugEigvf.append(values[7].strip())
            
            AugRate.append(values[9].strip())
    
    with open("../29/29-groundsatenergy.out", "r") as lines:
        header = lines.readline().strip() #header line
        lines.readline().strip() #header line
        lines.readline().strip() #header line
        #print(header.split("\t"))
        for i, line in enumerate(lines):
            values = line.strip().split()
            
            Configs.append(values)
    
    for conf in Configs:
        radShells = []
        
        for i, s in enumerate(AugShelli):
            if s not in radShells and AugShellf[i] == conf[1] and AugJJf[i] == conf[2] and AugEigvf[i] == conf[3]:
                radShells.append(s)
        
        RadTotRate = sum([float(rate) * (int(RadJJi[i]) + 1) for i, rate in enumerate(RadRate) if RadShelli[i] in radShells])
        AugTotRate = sum([float(rate) * (int(AugJJi[i]) + 1) for i, rate in enumerate(AugRate) if AugShelli[i] in radShells])
        
        AugSARate = sum([float(rate) * (int(AugJJi[i]) + 1) for i, rate in enumerate(AugRate) if AugShellf[i] == conf[1] and AugJJf[i] == conf[2] and AugEigvf[i] == conf[3]])
        
        if radShells:
            BR.append(AugSARate / (RadTotRate + AugTotRate))
        else:
            BR.append(0)
        
        print(conf[1] + " " + conf[2] + " " + conf[3] + ": BR = auger(" + str(AugSARate) + ") / radiative(" + str(RadTotRate) + ") + auger(" + str(AugTotRate) + ") = " + str(BR[-1]))
    
    for br in BR:
        print(br)
    
    
if __name__ == "__main__":
   main()