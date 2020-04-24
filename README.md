# **De Novo Genome Assembly Pipelines**
Re-usable block-builded,containerized pipelines for de novo genome assemblies.

<br>

## Welcome to a guide for building and maintaining the ready to go, containerized workflows of HCMR, for analyzing genome data.

![Banner](/banner.jpg)

<br>

> If you find yourself on this page, it means that either you are interested in building your own containerized images for genome analyses, and you have asked for access, either you are just curious in how the currently available ones were made. In both cases, before proceeding, make sure that you have checked what we have already build, understand how it works, and try to capture the purpose of this technology  firstly from the side of the user, [here](https://nellieangelova.github.io/DNGAW/). There is, of course, also the possibility of you being the one that has to update or maintain the images already build in HCMR. In all cases, keep reading and exploring!

<br>

## **About the containerized images**
> Let's start with some info and principles first. The pipeline analyses are made using [Snakemake](https://snakemake.readthedocs.io/en/stable/index.html) and containerized through [Singularity](https://sylabs.io/guides/3.5/user-guide/) with the help of the [Conda](https://docs.conda.io/en/latest/) package manager. If you are not familiar with any of these tools, this is the time for you to click on all three links above and get familiar with them, since they are the three pillars that make our technology work, and their understanding is fundamental for continuing. 

### **Why Snakemake?**
> Snakemake is a great workflow manager. The workflows created through it are somehow logic based: inputs and outputs of every step, with cool accessorizes such as parameters and wildcards that let you use multiple files at once and create directories or files or both with a single line of code. It's features make your pipeline's coding really easy to handle, quick and memoriable. Moreover, its special features such us creating directed acyclic graphs (DAGs) of the jobs composing a pipeline can save you a lot of time from working on graphs and images with other tools. It can even let you run a pipeline without actually running it, to detect problems, update files and search for new dependencies among the files, which by the way, understands by itself just by the names of the files. Cool ha? The workflows are portable, and can even work on clusters or the cloud. Importantly enough, it can be combined with the python package manager Conda, and the container platform Singularity.

### **Why Conda?**
> Conda is a package manager, created for the special task of creating isolated environments for important tools and jobs to be executed into, without interferring with the software of the host system. It actually creates mini-environments where you can run whatever you want, it installs the side packages and their dependencies you may need, switches between versions and installations according to your needs, and all that, without worrying you with compatability issues if handled right. It can create, save and switch between environments for many languages, such as Ruby and Java, but it was firstly and mainly created for Python users. And this can only mean, that its Python updates and features are the best in town. Now, if you have a Python based workflow manager that is composed by rules and tasks, and a Python based package manager that can create an isolated environment with everything this task may need to work, these two fellas can work together like chocolate and bannana. 


### **Why Singularity?**
> Singularity is a container platform, that lets you export the whole pipeline and its software in a single file, called an *image*. This image, can then be copied and executed in any environment that has Singularity installed and works with a linux based kernel. Everything that happens in the environment the image is spawning, affects the environment itself, meaning that the container image is using the host machine but does not interfer with its software part. Simply putted, the host does not need to have any of the tools, packages and softwares needed by the image to support it, which can only mean that your workflow is as portable as it can ever be. Moreover, since the whole pipeline is a single, autonomous image, it can also be submitted as a job in any cluster environment, which is and the actual point this project was build in the first place. The user of the image can run it in high performing environments and get his/her results as fast as possible, so he/she can continue researching in no time, and all that, without the need of dealing with compatability issues and various software installations that possibly require work by system administrators and changes in systems used by many users, such as in cluster environments. 


&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;![SCS](/scs.png)


<br>

## **About this project and its images**

<br>

> When it comes to genome analysis, the user can choose out of four different images, depending on his/her needs and of course, the available data. Here is a brief explanation for each one of them:

<br>

&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; 
![Images](/4-pipelines.png)

<br>

**A. QTQIllumina**

> QTQIllumina is an image suitable for those who want to run a simple quality control of their short reads data, before and after their trimming, when planning to use long and short reads combination for assembling a genome. It actually performs quality check, trimming of the raw data and quality check all over again. The results of the image are the following, located in the same directory your raw data live in:

1. One directory that contains the now trimmed short data, named *trimmomatic_output*
2. One directory that contains the results of the quality control before trimming the short reads, named *Multiq_raw_report*
3. One directory that contains the results of the quality control after trimming the short reads, named *Multiq_trimmed_report*

> Here are the tools used for this analysis:

| Tools       | Description        | Version |
| ------------- |:-------------:| -----:|
| fastqc     | Quality check |  |
| Multiqc     | Quality check |  |
| Trimmomatic    | Trimming     |   |

<br>

***
**B. QTQMinION**

> QTQMinION is an image suitable for those who want to run a simple quality control of their MinION data, before and after their trimming. It actually performs quality check, trimming of the raw data and quality check all over again. The results of the image are again, the following three, located in the same directory your MinION raw data live in:

1. One fastq file that contains the now trimmed data, named *porechop_output.fastq*
2. One directory that contains the results of the quality control before trimming, named *Nano_Raw_Report*
3. One directory that contains the results of the quality control after the trimming, named *Nano_Trimmed_Report*

> Here are the tools used for this analysis:

| Tools       | Description        | Version |
| ------------- |:-------------:| -----:|
| NanoPlot     | Quality check | 1.29.0 |
| Porechop    | Trimming     |  0.2.3 |

<br>

***


**C. LongGA**

> LongGA is an image suitable for those who want to run a whole genome assembly analysis from start to end, and got only long MinION reads at their disposal. Here are the main tools used by the pipeline, and all the results the image is spawning: 
>> QTQMinION is part of LongGA

1. One fastq file that contains the now trimmed data, named *porechop_output.fastq*
2. One directory that contains the results of the quality control before trimming, named *Nano_Raw_Report*
3. One directory that contains the results of the quality control after the trimming, named *Nano_Trimmed_Report*
4. One directory called *Assemblies*, which will actually contain all the assemblies made through the procedure and polishing steps alongside the final assembly. It contains the following subfolders and files:
* (DIR) Flye: A folder containing the new assembly called *assembly.fasta*, and other side files.
* (F) (x)_mapping.sam: A file generated through Minimap and used for polishing in Racon. x is a number and depends on the polishing rounds you've demanded in the configuration file. If for example your hyperparameter is 2, you will find two such files: 1_minimap.sam and 2_minimap.sam.
* (F) racon_(x).fasta: Polished assembly file generated through the Racon. (x) once again depends on the number of iterations. If the number of polishing rounds is 3, you will find 3 such files: racon_1.fasta, racon_2.fasta and racon_consensus.fasta. The last file is the final result of this part, that serves as input in Medaka.
* (F) racon_consensus.fasta.fai: A file for indexing of the latest racon fasta file.
* (F) racon_consensus.fasta.mmi: A file for indexing of the latest racon fasta file.
* (DIR) (Lineage_name): A folder containing information about the lineage used for the BUSCO analysis below. (e.g.*actinopterygii_odb9*)
* (F) The .tar form of the above lineage folder.
* (DIR) Medaka: A folder containing the final, fully polished assembly called *consensus.fasta*, and other side files.
5. One directory that contains the results of the BUSCO quality controls after each assembly creation, named *Busco_Results*. It contains two subfolders, *QA_1* and *QA_1*, for *assembly.fasta* and for *consensus.fasta* respectively.
6. One directory that contains the results of the QUAST quality controls after each assembly creation, named *Quast_Results*. It contains two subfolders, *QA_1* and *QA_1*, for *assembly.fasta* and for *consensus.fasta* respectively.
7. A pdf in the same directory with the image, with a directed acyclic graph (DAG) generated by it, which show the dependencies and input-outputs of the jobs of the image for you to check and use as you like.



| Tools       | Description        | Version |
| ------------- |:-------------:| -----:|
| NanoPlot     | Quality check | 1.29.0 |
| Porechop    | Trimming     | 0.2.3  |
| Flye   | Assembler  |  2.6 |
| Busco    | Quality check     |  3.0 (Internal Blast: v2.2) |
| Quast   | Quality check     |  5.0.2 |
| Racon   | Polishing     |  1.4.3 |
| Medaka   | Polishing    | 0.9.2  |


<br>

***

**D. Long_Short_GA**

> Long_Short_GA is an image suitable for those who want to run a whole genome assembly analysis from start to end, and got both long MinION reads and short Illumina reads at their disposal. Here are the tools and all the results the image is spawning: 
>> Long_GA and QTQIllumina are parts of Long_Short_GA
 

1. One directory that contains the now trimmed short data, named *trimmomatic_output*
2. A file called *KmerGenie.html* inside the trimmomatic_output generated by the KmerGenie program and informing about the estimated genome size and kmers of the data
3. One directory that contains the results of the quality control before trimming the short reads, named *Multiq_raw_report*
4. One directory that contains the results of the quality control after trimming the short reads, named *Multiq_trimmed_report*
5. One fastq file that contains the now trimmed data, named *porechop_output.fastq*
6. One directory that contains the results of the quality control before trimming, named *Nano_Raw_Report*
7. One directory that contains the results of the quality control after the trimming, named *Nano_Trimmed_Report*
8. One directory called *Assemblies*, which will actually contain all the assemblies made through the procedure and polishing steps alongside the final assembly. It contains the following subfolders and files:
* (DIR) Flye: A folder containing the new assembly called *assembly.fasta*, and other side files.
* (F) (x)_mapping.sam: A file generated through Minimap and used for polishing in Racon. x is a number and depends on the polishing rounds you've demanded in the configuration file. If for example your hyperparameter is 2, you will find two such files: 1_minimap.sam and 2_minimap.sam.
* (F) racon_(x).fasta: Polished assembly file generated through the Racon. (x) once again depends on the number of iterations. If the number of polishing rounds is 3, you will find 3 such files: racon_1.fasta, racon_2.fasta and racon_consensus.fasta. The last file is the final result of this part, that serves as input in Medaka.
* (F) racon_consensus.fasta.fai: A file for indexing of the latest racon fasta file.
* (F) racon_consensus.fasta.mmi: A file for indexing of the latest racon fasta file.
* (DIR) (Lineage_name): A folder containing information about the lineage used for the BUSCO analysis below. (e.g.*actinopterygii_odb9*)
* (F) The .tar form of the above lineage folder.
* (DIR) Medaka: A folder containing the final, fully polished assembly called *consensus.fasta*, and other side files.
9. One directory that contains the results of the BUSCO quality controls after each assembly creation, named *Busco_Results*. It contains two subfolders, *QA_1* and *QA_1*, for *assembly.fasta* and for *consensus.fasta* respectively.
10. One directory that contains the results of the QUAST quality controls after each assembly creation, named *Quast_Results*. It contains two subfolders, *QA_1* and *QA_1*, for *assembly.fasta* and for *consensus.fasta* respectively.
11. A pdf in the same directory with the image, with a directed acyclic graph (DAG) generated by it, which show the dependencies and input-outputs of the jobs of the image for you to check and use as you like.


| Tools       | Description        | Version |
| ------------- |:-------------:| -----:|
| fastqc     | Quality check |  |
| Trimmomatic     | Trimming |  |
| Multiqc     | Quality check |  |
| NanoPlot     | Quality check | 1.29.0 |
| Porechop    | Trimming     | 0.2.3  |
| Flye   | Assembler  |  2.6 (Internal: KmerGenie: v1.70.16)|
| Busco    | Quality check     |  3.0 (Internal Blast: v2.2) |
| Quast   | Quality check     |  5.0.2 |
| Racon   | Polishing     |  1.4.3 |
| Medaka   | Polishing    | 0.9.2  |
| Pilon  | Polishing    |  1.23 (Internal: Samtools: v1.7, Minimap2: v2.17)|

<br>

***


## **A LEGO logic**

The main idea behind the composition of this project, was to build something out of individual pieces that can be used solely, maintained easily, and be generally *independent*. And here comes Conda.
<br>
[LEGO](/lego.jpg)

```
singularity run <image.simg>
```

***
#### Please credit accordingly:
> Author: Nellie Angelova, Bioinformatician, Hellenic Centre for Marine Research (HCMR) <br>
> Coordinator/ Supervisor: Tereza Manoussaki, Researcher, Hellenic Centre for Marine Research (HCMR)

