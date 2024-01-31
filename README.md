# ASFV_Pipeline
## Dependencies
Pipeline has been created using a docker image. Singularity may also be used. Please install one of these programs.

1) Docker
https://docs.docker.com/engine/install/

	OR

2) Singularity
https://docs.sylabs.io/guides/3.0/user-guide/installation.html

## Installation
Transfer the designated folders and files from Github to the directory where you intend to execute the pipeline.
### Files
1) MetaDataStorage.csv 
2) MetadataNew.csv
3) Graphic.txt
4) BiotypingGenomes.csv
5) LinuxPipeLine_v01_container.py
6) p72GenotypingData.csv

### Folder
./BlastOutput/

### Download the Docker Image
The Image can be downloaded and formatted for docker or singuilarity using the following commands:

#### Docker
First, open Docker Desktop. Then in Powershell enter the following command:
```
 docker pull garadock/asfv_match_assembly_pipeline:v02
 ```
#### Singularity
Open Powershell and enter the following command:
```
singularity pull ASFV_match_assembly_pipeline_v02.sif docker://garadock/asfv_match_assembly_pipeline:v02
```

## Instructions for use

### Set up MetadataNew.csv
1) In MetadataNew.csv, fill out the followings columns (L1R1,L1R2,L2R1,L2R2,L3R1,L3R2,L4R1,L4R2)

    a) These columns refer to the Lane "L" and Read "R" number.
   
    b) Please put the path to the fastq.gz file. Example:
   
        L1R1 = Input/SampleAA_R1_001.fastq.gz
   
        L1R2 = Input/SampleAA_R2_001.fastq.gz
   
    c) L1R1 and L1R2 must have reads, the remaining columns can remain blank if the data files were not generated.
   
3) In MetadataNew.csv, fill out the followings columns (Email,collection_date,country,location,host,tissue_type,collected_by,isolate)
   
    a) These columns are less important, they can contain any string.
   
5) Keep column "Project_ID (internal)" blank

### Running the Script (Docker)
Enter the directory that contains the fastq.gz files and the files and folders detailed in part 1 of the installation section.

#### Windows
1) Open docker desktop

2) Open powershell and paste the following commands:

    ```powershell
    $LinHomeMount = $HOME.substring(2).replace('\','/')
    $LinPWD = $pwd.path.substring(2).replace('\','/')
    docker run -it -v ${HOME}:$LinHomeMount -w $LinPWD garadock/asfv_match_assembly_pipeline:v02 /bin/bash
    ```

   ```
   python3 LinuxPipeLine_v01_container.py
   ```

#### Mac
Coming Soon
#### Linux
Coming Soon

### Running the Script (Singularity)
Example .sh script:
```
#!/usr/bin/env bash
#SBATCH --ntasks=#
#SBATCH --nodes=#
#SBATCH --partition=string
#SBATCH --cpus-per-task=## #Comment: On Line 23 of the LinuxPipeLine_v01_container.py, the number of threads is set to 10, this can be increased to the number of CPUS indicated here on the cluster.
#SBATCH --job-name=jobname
#SBATCH --time=##
#SBATCH --output=log.log


module load singularity


File_DIR=/your/directory/conatining/files/

export SINGULARITY_BIND=$File_DIR

singularity exec --writable-tmpfs "/path/to/the/container/ASFV_match_assembly_pipeline_v02.sif" \
python3 LinuxPipeLine_v01_container.py
```
