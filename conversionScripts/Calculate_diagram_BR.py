import sys

RadShelli = []
RadJJi = []
RadShellf = []
RadJJf = []
RadEigvf = []
RadRate = []

AugShelli = []
AugJJi = []
AugShellf = []
AugJJf = []
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
            RadEigvf.append(values[7].strip())
            
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
            
            AugShellf.append(values[5].strip())
            AugJJf.append(values[6].strip())
            
            AugRate.append(values[9].strip())
    
    with open("../29/29-grounddiagenergy.out", "r") as lines:
        header = lines.readline().strip() #header line
        lines.readline().strip() #header line
        lines.readline().strip() #header line
        #print(header.split("\t"))
        for i, line in enumerate(lines):
            values = line.strip().split()
            
            Configs.append(values)
    
    for conf in Configs:
        radShells = []
        
        for i, s in enumerate(RadShelli):
            if s not in radShells and RadShellf[i] == conf[1] and RadJJf[i] == conf[2] and RadEigvf[i] == conf[3]:
                radShells.append(s)
        
        RadTotRate = sum([float(rate) * (int(RadJJi[i]) + 1) for i, rate in enumerate(RadRate) if RadShelli[i] in radShells])
        AugTotRate = sum([float(rate) * (int(AugJJi[i]) + 1) for i, rate in enumerate(AugRate) if AugShelli[i] in radShells])
        
        RadSARate = sum([float(rate) * (int(RadJJi[i]) + 1) for i, rate in enumerate(RadRate) if RadShellf[i] == conf[1] and RadJJf[i] == conf[2] and RadEigvf[i] == conf[3]])
        
        if radShells:
            BR.append(RadSARate / (RadTotRate + AugTotRate))
        else:
            BR.append(0)
        
        print(conf[1] + " " + conf[2] + " " + conf[3] + ": BR = radiative(" + str(RadSARate) + ") / radiative(" + str(RadTotRate) + ") + auger(" + str(AugTotRate) + ") = " + str(BR[-1]))
    
    for br in BR:
        print(br)
    
    
if __name__ == "__main__":
   main()