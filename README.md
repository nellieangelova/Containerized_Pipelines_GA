# **Containerized Pipelines for Genome Assembly Analyses**
Re-usable, block-builded, containerized pipelines for de-novo genome assembly construction analyses.

<br>

## Welcome to a guide for building and maintaining the ready to go, containerized workflows of HCMR, for analyzing genome data.

![Banner](/banner.jpg)

<br>

> If you find yourself on this page, it means that either you are interested in building your own containerized images for genome analyses, and you have asked for access, either you are just curious in how the currently available ones were made. In both cases, before proceeding, make sure that you have checked what we have already build, understand how it works, and try to capture the purpose of this technology  firstly from the side of the user, [here](https://nellieangelova.github.io/DNGAW/). There is, of course, also the possibility of you being the one that has to update or maintain the images already build in HCMR. In all cases, keep reading and exploring!

<br>

## **About the containerized images**
> Let's start with some info and principles first. The pipeline analyses are made using [Snakemake](https://snakemake.readthedocs.io/en/stable/index.html) and containerized through [Singularity](https://sylabs.io/guides/3.5/user-guide/) with the help of the [Conda](https://docs.conda.io/en/latest/) package manager. If you are not familiar with any of these tools, this is the time for you to click on all three links above and get familiar with them, since they are the three pillars that make this technology work, and their understanding is fundamental for continuing. 

### **Why Snakemake?**
> Snakemake is a great workflow manager. The workflows created through it are somehow logic based: inputs and outputs of every step, with accessorizes such as parameters and wildcards, that let you use multiple files at once and create directories (or files, or both) with a single line of code. Its features make your pipeline's coding really easy to handle, quick and memoriable. Moreover, its special features such us creating directed acyclic graphs (DAGs) of the jobs composing a pipeline can save you a lot of time from working on graphs and images with other tools. It can even let you run a pipeline without actually running it, to detect problems, update files and search for new dependencies among the files, which by the way, understands by itself just by the their names. The workflows are portable, and can even work on clusters or the cloud. Importantly enough, it can be combined with the python package manager Conda, and the container platform Singularity.

### **Why Conda?**
> Conda is a package manager, created for the special task of creating isolated environments for important tools and jobs to be executed into, without interferring with the software of the host system. It actually creates mini-environments where you can run whatever you want, it installs the side packages and their dependencies you may need, switches between versions and installations according to your needs, and all that, without worrying you with compatability issues if handled right. It can create, save and switch between environments for many languages, such as Ruby and Java, but it was firstly and mainly created for Python users. And this can only mean, that its Python updates and features are carefully maintained. Now, if you have a Python based workflow manager that is composed by rules and tasks, and a Python based package manager that can create an isolated environment with everything this task may need to work, you've got yourself some powerful technology to work on. 


### **Why Singularity?**
> Singularity is a container platform, that lets you export the whole pipeline and its software in a single file, called an *image*. This image, can then be copied and executed in any environment that has Singularity installed and works with a linux based kernel. Everything that happens in the environment the image is spawning, affects the environment itself, meaning that the container image is using the host machine but does not interfer with its software part. Simply putted, the host does not need to have any of the tools, packages and softwares needed by the image to support it, which can only mean that your workflow is as portable as it can ever be. Moreover, since the whole pipeline is a single, autonomous image, it can also be submitted as a job in any cluster environment, which is and the actual point this project was build in the first place. The user of the image can run it in high performing environments and get his/her results as fast as possible, so he/she can continue researching in no time, and all that, without the need of dealing with compatability issues and various software installations that possibly require work by system administrators and changes in systems used by many users, such as in cluster environments. 


&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;![SCS](/scs.png)


<br>

## **About this project and its images**
<br>

> When it comes to genome analysis, the user can choose out of four different images, depending on his/her needs and of course, the available data. If you've checked the userguide, you've probably already come across this image below. As you understand, the available analyses are parts of a bigger analysis and the whole project is somehow block-based.

<br>

&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; 
![Images](/4-pipelines.png)

<br>

## **A LEGO logic**

&nbsp; &nbsp; &nbsp; &nbsp; &nbsp;&nbsp; &nbsp; &nbsp; &nbsp; &nbsp;&nbsp; &nbsp; &nbsp; &nbsp; &nbsp;&nbsp; &nbsp; &nbsp; &nbsp; &nbsp;![LEGO](/lego.jpg)

> The main idea behind the composition of this project, was to build something out of individual pieces that could be used solely, maintained easily, and be generally *independent*. And here comes Conda to the rescue. An analysis such as a genome assembly analysis is composed out of many steps, which can be seen as independent jobs. These jobs are the blocks, that all together create the main workflow. Each job has an input and gives an output, from another job or from the user, to another job and/or the user respectively. In Snakemake terms, these jobs can be refered as *rules*. 
> By using Singularity, a universal conda environment is created, which can host independent mini-environments and "pull their strings". With the help of conda, each rule has its own environment to run into. This environment includes a specific python version, if needed by the software it uses, and the softwares themselves, which are the basic tools needed by the rule to generate the desired output. For example, in the case of the Polishing step in the genome assembly analysis, three tools are needed to polish a Flye made assembly: Minimap2, Racon and Medaka. Each environment uses by default the version of Python the base conda environment was using during their creation. In this case, this is Python 3.6.10. If for example Racon was running only for versions of Python <3.5, we could and probably would make the Polishing environment to turn to that release so the rule can be executed. This environment is opened when the rule is about to be executed, and closes when the output is formed, without the universal environment being affected by any of these changes. 

> The isolation of each process that can be seen as a mini-workflow, in an environment with specific tunnings, is very helpful when one needs to organize his/her workflow. It is very difficult for a developer to track all the changes made in a big, one-piece project with lots of data, packages and code. The **Divide and Conquer** logic here, where a main task is splitted in many sub-tasks running independently, gives a developer the opportunity to:
* Organize better his/her code depending on the chunk, so having in fact better control over the project.
* Re-use pieces of code and mini-workflows, that can easily stand on their own as images for specific tasks or participate in the creation of future images.
* Maintain his/her project much easier, since updating or changing the release of one program won't mean changing the whole project itself. Just a piece of it, which can be downloaded, changed, and pushed back in its place in no time. Like a piece of Jenga!

<br>

&nbsp; &nbsp; &nbsp; &nbsp; &nbsp;&nbsp; &nbsp; &nbsp; &nbsp; &nbsp;&nbsp; &nbsp; &nbsp; &nbsp; &nbsp;&nbsp; &nbsp; &nbsp; &nbsp; &nbsp;&nbsp; &nbsp; &nbsp; &nbsp; &nbsp;&nbsp; &nbsp; &nbsp; &nbsp; &nbsp;&nbsp; &nbsp; &nbsp; &nbsp; &nbsp;&nbsp; &nbsp; &nbsp; &nbsp; &nbsp;![JENGA](/jenga.jpg)

## **Diving Deeper into the creation and composition of each environment and subsequently the images**
<br>

> If you are not familliar with the basic principles of Python, Conda, Snakemake and Singularity yet, this is the right time for you to **worry**. 
>> Each script and code mentioned here can be found in the folders and subfolders of the project in the current GitHub repository, fully commented.

### **Step 1: Creating a general environment of Singularity, and subsequently an image**
>In order for the code to live somewhere, and this somewhere to be portable, you need Singularity. Now, according to its documentation, you can build an environment that can be a whole (almost) operating system. Here,  a **conda super-environment** is created. This new environment, can work wherever you run the final Singularity image we are going to build. The Singularity file is called *"<name_of_the_image>Singularity"* and it is written in a simple text form, as Singularity has its own "language" for expressing things. When you open such a file, you will see some specific blocks of code. 
* **Header**: The header ("Bootstrap", "From") is at the top of the file, and defines the base you want to use to build the container on to. It briefly describes the os or environemnt you wish to build everything in to.
* **%files**: The files section uses the traditional copy (cp) command. Anything declared here, from files to folders, is included **in** the final image. So every single script or anything else may be necessary for the project to run, but we don't want the user to have anything to do with it, is included in the container a priori. This files can be found at the **root (/)** of the image. 
* **%environment**: In this section, variables needed during the *runtime* are declared. Since no one has interaction with the image during runtime, things such as adding something to the $PATH a priori is essential.
* **%post**: The important thing to note here is that a Singularity image, may work like an independent system, but uses the Linux-based kernel of it's host. Since a user inside the image cannot and should not have sudo privileges on the host system, some basic commands such as activating environments should be written in the *bashrc* of the container during *build* time, while one builds it's images as an admin on his own machine. 
* **%runscript**: The runscript section includes everything you want to run during runtime. When the image is ready, and run via the *singularity run* command, everything included here will be executed  progressively. In the case of this project, the user is asked to just run the image as a job into his cluster, with the simple command:

```
singularity run <image.simg>
```

>When he does, he actually activates some code from the inside of the image. For example, if he runs the LGA workflow, he activates the following from the *LGASingularity* file:
```
source activate /opt/conda/envs/LQTQ
snakemake -j --snakefile /LGASnakefile --use-conda --nolock --quiet --keep-going
snakemake -j --snakefile /LGASnakefile --dag | dot -Tpdf > DAG.pdf
snakemake -j --snakefile /LGASnakefile --summary > Summary.txt

rm -R config
```
>The first row activates the first sub-environment, which also has snakemake installed, to begin the process. The second one is actually the command that initializes the workflow.
>> Braking down the snakemake command used:
* j: The number of available cores. If the number is ommited, it is determined by Snakemake as the number of available CPU cores the machines has (for parallelism, discussed later).
* snakefile: The file which includes the rule definitions. Included during the %files stage while building the image.
* use-conda: Make conda use the specific environment mentioned in each rule for its execution.
* nolock: Does not lock the working directory if something goes wrong, so the input files can be accesable for a future re-run.
* quiet: Keeping Snakemake from outputting excessive information about each rule's execution.
* keep-going: Go on with independent jobs if a job fails.
> The third row, makes snakemake generate a .pdf with the dependencies between the tasks, as an informative diagram for the user.
> The forth row creates a summary with the files generated by the pipeline, so the user can check whether all the outputs are created as expected or errors have occured. 
> The fifth row deletes a temporary folder from the user's workdir.

<br>

> To build the image, go to the directory your files and your Singularity definition file lives in, and just run:

```
sudo singularity build <name_of_the_image>.simg <name_of_the_Singularity_Definition_File>

```
> This will result in a .simg file in the same directory, which can be transfered and copied anywhere you wish to run it.

### **Step 2: The Snakemake Pipeline**
> But what about the pipeline and its rules?
> The pipeline needs also its definition file, mentioned above, here called *"<name_of_the_image>Snakefile"*. It is copied inside the image and it lives in its root, so it can be called from there during runtime. This file holds the definition of each of the rules. Each rule has input(s), output(s), an environment to use, and a script with the code that needs to be executed. It can also have some parameters, and a number of threads to use. Each environment used by each rule has its own definition file in YAML format, which is also included in the image during the %files stage. The Snakefile also contains the declaration of an existing configuration file to use, *"config.yaml"*, which contains the hyperparameters the workflow is going to use. This file is given to the user and has to be filled by him according to his needs, so it does not live inside the image, but should be copied by him in the same directory he copies the image into, for snakemake to be able to find it. Additionally, it has a rule called *all*, which has as inputs all the outputs the image has to give to the user, and assures that everything is generated smoothly during runtime. If an output file fails, the rule *all* fails, and thus the image.

> **Multi-threading**
> The *threads* field is one of the most important ones, since it provides the needed information for Snakemake to schedule the jobs it submits to the available cores of the host machine. The j parameter of Snakemake, in the Singularity definition file, if ommited, lets Snakemake use all the available cores of the machine. The threads that each rule needs are be a portion of these available cores. The number is considered as the *maximum* cores a rule needs, and thus it is reduced depending on the available cores of the system. If the threads field is ommited, Snakemake uses its default number of 1 thread per rule, and submits all the rules together, overloading the environment. Deciding the number of threads a rule may need for optimizing a pipeline, takes time and needs tunning. 


<br>

> One who needs to maintain or update an environment, has to do the following:
1. Download the .yml file of the environment of interest.
2. Create an environment from it on his/her machine.
3. Make the changes.
4. Export the environment in a new .yml file, which will take the place of the previous one.
5. Build again the affected images with the new .yml file.
6. Run the images to assure that everything works smoothly before giving them to the community.
7. Update the webpages and related information lists to keep track of the changes and inform the community about them.

## **The environments and their use**
> Here is a detailed list of all currently available environments and the images they are hosted into.
>> If you took the time to check the available images in the userguide, as suggested, you already know that there are four of them (SQTQ, LQTQ, LGA and LSGA), and each one is more or less a part of or an extension of another. 

<br>

**SQTQ**:
> SQTQ serves in two images (SQTQ, and LSGA) and five rules (fastq_c, multic_c, trimming_s, fastq_c_t, multiq_c_t). It contains programs that check the quality of short raw Illumina data, trim them, and then check the quality of the trimming.

| Tools       | Description        | Version |
| ------------- |:-------------:| -----:|
| Python   | Version of the python version of the base Conda used, during build time. | 3.6.10 |
| Snakemake     | Since this environment serves and as a standalone for an image, it has to serve as base and include Snakemake. SQTQ is activated for SQTQ and LSGA images.| 5.3.0|
| fastqc     | Performs quality check of the Illumina data.| 0.11.8| 
| Multiqc     | Creates a summary of the fastqc runs for each Illumina file. | 1.6 |
| Trimmomatic    | Trimms the Illumina raw data.     |  0.39 |

<br>

**LQTQ**:
> LQTQ serves in three images (LQTQ, LGA and LSGA) and three rules (nanoq_c, trimming_l, nanoq_c_t). It contains programs that check the quality of long raw MinION data, trim them, and then check the quality of the trimming.


| Tools       | Description        | Version |
| ------------- |:-------------:| -----:|
| Python   | Version of the python version of the base Conda used, during build time. | 3.6.10 |
| Snakemake     | Since this environment serves and as a standalone for an image, it has to serve as base and include Snakemake. It serves as base also for the rest two images. LQTQ is activated for LQTQ and LGA images.| 5.3.0 |
| NanoPlot     | Creates a summary of the quality check run for the MinION file. | 1.29.0 |
| Porechop    | Trims the MinION data.  |  0.2.3 |

<br>

**FlyeAss**:
> FlyeAss serves in two images (LGA and LSGA) and one rule (FlyeAssGenie). It contains programs that estimate the genome size with the help of short Illumina data, if any, and create the first assembly out of the MinION data.

|Tools       | Description        | Version |
| ------------- |:-------------:| -----:|
|Python   | An older version needed by KmerGenie. | 2.7.17 |
|Flye   | MinION Assembler |  2.6 |
|KmerGenie  | Genome size estimation for the upcoming assembly, out of Illumina data (if any). |1.70.16|

<br>

**QAssembly**:
> QAssembly serves in two images (LGA and LSGA) and three rules (QA_1,QA_2). It contains programs that estimate the quality of a given assembly.

|Tools       | Description        | Version |
| ------------- |:-------------:| -----:|
| Python   | Version of the python version of the base Conda used, during build time. | 3.6.10 |
| Busco    | Quality check of a given assembly.     |  3.0.1 (Internal Blast: v2.2, due to multithreading errors of newer versions.) |
| Quast   | Quality check of a given assembly.    |  5.0.2 |

<br>

**Polishing**:
> Polishing serves in two images (LGA and LSGA) and one rule (Polishing). It contains programs that polish a given assembly file in rounds.

|Tools       | Description        | Version |
| ------------- |:-------------:| -----:|
| Python  | Version of the python version of the base Conda used, during build time. | 3.6.10 |
| Racon   | Polishes a given assembly file.   |  1.4.3 |
| Medaka  | Creates some files (e.g. bam files) needed by Racon for polishing.    | 0.9.2  |

<br>

**Piloning**:
> Piloning serves in just one image (LSGA) and one rule (Piloning). It contains programs that polish a given assembly file through correcting errors by using short Illumina reads.

| Tools       | Description        | Version |
| ------------- |:-------------:| -----:|
| Python  | Version of the python version of the base Conda used, during build time. | 3.6.10 |
| Pilon  | Polishing a given assembly file with both MinION and Illumina reads for creating a final consensus assembly.   |1.23 (Internal: Samtools: v1.9, Minimap2: v2.17)|
| Busco    | Quality check of a given assembly.     |  3.0.1 (Internal Blast: v2.2, due to multithreading errors of newer versions.) |
| Quast   | Quality check of a given assembly.    |  5.0.2 |


<br>

> To understand how each workflow works, you should check its Singularity definition file first, to see what it contains and the parameters of Snakemake. Then, check its Snakemake definition file, to understand the flow of the inputs and outputs for generating the desired output(s), and also the threading. Open the script responsible for each rule, and check the tools of its environment. 
>> When developing something, the first and last thing you should do is walking in the shoes of your future users. Create the configuration files they will need, inform them about the additional steps they should make before and after running your code, and keep your documentaries and subsequently your users updated. 

***
#### Please credit accordingly:
> Author: Nellie Angelova, Bioinformatician, Hellenic Centre for Marine Research (HCMR) <br>
> Coordinator/ Supervisor: Tereza Manoussaki, Researcher, Hellenic Centre for Marine Research (HCMR) <br>
> Contact: n.angelova@hcmr.gr

