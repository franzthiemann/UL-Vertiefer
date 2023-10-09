#!/usr/bin/env python3
import subprocess
import re
import os
import pandas as pd
vmdfile = "Path to vmd template file for rendering"
vmdrc = "Path to vmdrc configuration file"
csvfile = "Path to the result summary file"

class pEDAproperties:
    def __init__(self, runfile, folder, type, NOCVs, structure, full, surface, mol, unrestricted) -> None:
        self.runfile = runfile
        self.folder = folder
        self.type = type
        self.NOCVs = NOCVs
        self.structure = structure
        self.full = full
        self.surface = surface
        self.mol = mol
        self.unrestricted = unrestricted
    
    def __str__(self):
         out = f"Runfile: {self.runfile} --> Type: {self.type}\n"
         out = out + f"Results in Folder: {self.folder}\n"
         if self.NOCVs > 0: out = out + f"{str(self.NOCVs)} NOCVs generated"
         out = out + f"Linked Vasp Files: {self.full}; {self.surface}; {self.mol}"
         return out


def analyzeRunfile(file):
    print("Found Run File " +  file)
    f = open(file, "r")
    lines = f.readlines()
    for line in lines:
        if re.match(r"^pEDA\s", line):
                # Command line
                folder = re.search(r"\s-f\s(\S+)\s", line)[1]
                type = re.search(r"\s--task\s(\S+)\s", line)[1]
                if type == "NOCV" and re.search(r"\s--NOCVs\s(\S+)\s", line):
                     NOCVs = int(re.search(r"\s--NOCVs\s(\S+)\s", line)[1])
                else:
                     NOCVs = 0
                structure = re.search(r"(?:^|\s)[^\s-]+\s([^\s-]+)(?=\s|$)", line)[1]
                full = re.search(r"--vasp_outfile\s(\S+)\s", line)[1]
                surface = re.search(r"--relaxed_surf_vasp\s(\S+)\s", line)[1]
                mol = re.search(r"--relaxed_mol_vasp\s(\S+)\s", line)[1]
                
                if re.search(r"-u True", line):
                     unrestricted = True
                else:
                     unrestricted = False
                properties = pEDAproperties(file, folder, type, NOCVs, structure, full, surface, mol, unrestricted)
                return properties
    return False

def printFolderName(properties):
     print(f"*******{properties.folder}*******")

def evalpEDA(properties):
     # pEDA_eval --surf_ams $i/surface --mol_ams $i/mol --vasp_outfile adsorbat.vasp-out --relaxed_surf_vasp surface.vasp-out 
     # --relaxed_mol_vasp mol.vasp-out -u False --NOCVs 3 $i
     cmdarray = ["pEDA_eval "]
     if properties.unrestricted:
          cmdarray.append("--unrestricted True ")
     else:
         cmdarray.append("--unrestricted False ")
     cmdarray[1] += f"--surf_ams {properties.folder}/surface "
     cmdarray[1] += f"--mol_ams {properties.folder}/mol "
     cmdarray[1] += f"--vasp_outfile {properties.full} "
     cmdarray[1] += f"--relaxed_surf_vasp {properties.surface} "
     cmdarray[1] += f"--relaxed_mol_vasp {properties.mol} "
     if properties.NOCVs > 0:
           cmdarray[1] += f"--NOCVs {properties.NOCVs} "
         # cmdarray.append(f"--NOCVs {properties.NOCVs}")
     cmdarray[1] += properties.folder
     # print(cmdarray)
     command = cmdarray[0] + cmdarray [1]
     print(command)
     output = os.system(command + ">/dev/null")
     # print(output)
     # load csv
     processCSV(properties)

def processCSV(properties):
     data = pd.read_csv(properties.folder + ".csv")
     transposedData = data.transpose()
     # check for alpha / beta mentions in the date
     columns = transposedData.iloc[0, :].tolist()
     for nr, column in enumerate(columns):
          # print(column)
          if re.search(r"(_alpha_|_beta_|->)", column):
               transposedData.drop([nr], axis=1)
               # print(column)
     absf = os.path.abspath(properties.folder)
     transposedData.insert(0, "file", [absf, absf, absf])
     if os.path.isfile(csvfile):
          with open(csvfile, "a") as file:
               file.write(transposedData.to_csv(index=False, header=False))
     else:
          transposedData.to_csv(csvfile, index=False)


     # print(data)

def rendercube(properties):
     # get cube files
     files = [f for f in os.listdir(properties.folder) if os.path.isfile(properties.folder + "/" + f)]
     #create render scripts
     for file in files:
          if re.match(r".+\.cube$", file):
               renderscript = open(properties.folder + "/" + file + ".vmd", "w+")
               renderscript_str = "mol new " + properties.folder + "/" + file + " type cube first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all\n"
               renderscript_str += "play " + vmdrc + "\n"
               renderscript_str += "play " + vmdfile + "\n"
               renderscript_str += "render X3D " + properties.folder + "/" + file + ".x3d\n"
               renderscript_str += "quit"
               renderscript.write(renderscript_str)
               renderscript.close()
               # render the files
               command = "vmd -e " + properties.folder + "/" + file + ".vmd"
               os.system(command + ">/dev/null")
               print("rendered file " + file)


               

def folderwalker():
    files = [f for f in os.listdir('.') if os.path.isfile(f)]
    for file in files:
        if re.match(r".+\.run", file):
            properties = analyzeRunfile(file)
            # print(properties)
            printFolderName(properties)
            pEDARes = evalpEDA(properties)
            if properties.type == "NOCV" and properties.NOCVs > 0:
               rendercube(properties)

folderwalker()
