
import subprocess
import os

def quals(datadir,datanew, t):

    ''' This function takes raw data and performs quality control.
        datadir must be the directory that contains the data'''

    if t:
        os.system("NanoPlot --fastq " + snakemake.input[0] + " -o " + datanew)

    else:
        fqf = subprocess.run(["ls","-R", datadir], stdout=subprocess.PIPE, check = True).stdout
        lof = [x.strip() for x in str(fqf).split('\\n')]
        for f in lof:
            if ("fastq") in f:
                os.system("NanoPlot --fastq " +datadir+"/"+f + " -o "+ datanew)

    return


quals(snakemake.input[0],snakemake.output[0], snakemake.params["Trimmed"])
