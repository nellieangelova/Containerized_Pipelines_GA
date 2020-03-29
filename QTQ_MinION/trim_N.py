
import os


def trimming(datadir):

    ''' This function takes raw data and performs trimming.
        datadir must be the directory that contains the data.
        The output files are saved in a subfolder called "porechop_output".'''

    for f in os.listdir(datadir):
        if ("fastq" in f):
            os.system("porechop -i " + datadir +"/"+f +  " -o " + datadir+"/porechop_output.fastq -v 0  --discard_middle")

    return


trimming(snakemake.input[0])