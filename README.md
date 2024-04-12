# ASFV Pipeline
## Dependencies
Pipeline has been created using a docker image. Singularity may also be used. Please install one of these programs.

1) Docker Desktop
https://docs.docker.com/engine/install/

	OR

2) Singularity
https://docs.sylabs.io/guides/3.0/user-guide/installation.html

## Installation
Transfer the designated folders and files from Github to the directory where you intend to execute the pipeline. *It is recommended that you set this up in a new folder with files backed up.*
### Files
1) MetaDataStorage.csv 
2) MetadataNew.csv
3) ASFV_Pipeline.py
4) DeNovo_Functions.py

### Download the Docker Image
The Image can be downloaded and formatted for docker or singuilarity using the following commands:

#### Docker
First, open Docker Desktop. Then in Powershell enter the following command:
```powershell
docker pull garadock/asfv_denovo_assembly_pipeline:v05
```
#### Singularity
Open the terminal and enter the following command:
```terminal
singularity build --fakeroot ASFV_denovo_assembly_pipeline_v05.sif docker://garadock/asfv_denovo_assembly_pipeline:v05
```

## Instructions for use

### Set up MetadataNew.csv
MetadataNew.csv can be set up to run a single run or multiple consecutive runs. For each sample, fill in the row based on the following instructions:

1) In MetadataNew.csv, fill out the followings columns for Illumina reads (L1R1,L1R2,L2R1,L2R2,L3R1,L3R2,L4R1,L4R2)

    a) These columns refer to the Lane "L" and Read "R" number.

    b) Please put the paths to the fastq.gz files. Example:

        L1R1 = Input/SampleAA_R1_001.fastq.gz

        L1R2 = Input/SampleAA_R2_001.fastq.gz

    c) L1R1 and L1R2 must have reads, the remaining columns can remain blank if the data files were not generated.
    
2) (Optional) In MetadataNew.csv, fill out the following column for Nanopore reads (MinionDirectory).

    a) This column refers to the directory with reads generated from Oxford Nanopore Technologies. If no nanopore reads were generated, leave blank.

    b) Please put the path to the folder containing the reads. Example:

        Input/Minion/

    c) Be sure in to include the last forward slash "/"

3) In MetadataNew.csv, fill out the followings columns (Email,collection_date,country,location,host,tissue_type,collected_by,isolate)
   
    a) These columns are less important, they can contain any string.
   
4) (Optional) In MetadataNew.csv, you can provide a custom project name. The output will include date/time and a short custom string for identification purposes within the program. The input will be adjusted to be computer friendly. If one is not provided a randomly generated string and date will be used as a project name.

#### Example of Filling Out MetaDataNew.csv
This is an example of a properly filled out MetaDataNew.csv with two rows. You can click on the image to zoom in. Notice that forward slashes "/" are used instead of backslashes "\".
![image](https://github.com/Global-ASFV-Research-Alliance/ASFV_Pipeline/assets/103464896/247b41a9-e4d5-4795-972f-4172f8d9c991)

### Running the Script (Docker)
Enter the directory that contains the path to fastq.gz files, Minion directory (optional), and the files and folders detailed in part 1 of the installation section. Example:

```powershell
cd Documents
cd ASFV_Pipeline
```

#### Windows
1) Open docker desktop

2) Open powershell and paste the following commands:

    ```powershell
    $LinHomeMount = $HOME.substring(2).replace('\','/')
    $LinPWD = $pwd.path.substring(2).replace('\','/')
    docker run -it -v ${HOME}:$LinHomeMount -w $LinPWD garadock/asfv_denovo_assembly_pipeline:v05 /bin/bash
    ```

    ```powershell
    python3 ASFV_Pipeline.py
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
#SBATCH --cpus-per-task=## #Comment: On Line 23 of the DeNovo05.py, the number of threads is set to 10, this can be increased to the number of CPUS indicated here on the cluster.
#SBATCH --job-name=CustomJobName
#SBATCH --time=##
#SBATCH --output=CustomLogName.log


module load singularity


File_DIR=/your/directory/containing/files/

export SINGULARITY_BIND=$File_DIR

singularity exec "/path/to/the/container/ASFV_denovo_assembly_pipeline_v05.sif" \
python3 ASFV_Pipeline.py
```
## Troubleshooting, Bugs, and Known Problems
### My Program Dies Early On
1)	Make sure MetaDataNew.csv is set up properly. This includes changing all "\" to "/".

## Changelog
### ReadMe Update - Version 0.5a
No update of code required, this update solely updates the ReadMe.md file to improve clarity of certain subjects.
#### Minor Updates
1)	Example table included to clarify how MetaDataNew.csv should be constructed.

### Version 0.5 - DeNovo Release
#### Major Updates
1)	PLEASE RUN IN NEW FOLDER WITH BACKUP COPIES OF FILES ELSEWHERE. A new script for deleting temporary files has been introduced – this should help save space, but we would like to make sure the risk of deleting files unintentionally is minimized.
2)	De novo pipeline introduced – requires Illumina and Nanopore reads – greater accuracy provided along with no match needed. The old pipeline will be referred to as the “Illumina Pipeline”. Minion directory can be specified in new MetaDataNew.csv column.
3)	Initial mapping to Swine included to get rid of potential contamination. Gathered reads are those that do not map to Swine and those that map to ASFVG. Indexed swine genome is included in docker image.
4)	Improvements to the Identification of repeat regions.

#### Minor Updates
1)	Multiple runs can be included on a single MetaDataNew.csv file – with each row for a different run.
2)	Check for duplicate/absent lanes introduced – specific error will be indicated in the MetaDataStorage.csv file.
3)	Custom project names can be included in the MetaDataNew.csv file (Just type one in!). It will automatically be adjusted to be computer friendly.
4)	Some identified errors are now printed out in the MetaDataStorage.csv, the run will be stopped, and the next one will continue. Temp files are not deleted to allow for analysis.
5)	Updated explanation of pipeline in Summary readout.
6)	Fewer files needed in directory to run the process.

#### Illumina Only Updates
1)	Assembled scaffolds that show high % identity to swine are removed.
2)	For finding best match, only contigs greater than or equal to 1250 in length will be considered.
3)	De Novo Scaffolding mapped to predicted reference included in outputs.

## References
1) Bankevich, A. et al. SPAdes: a new genome assembly algorithm and its applications to single-cell sequencing. J Comput Biol. 2012 May;19(5):455-77. doi: 10.1089/cmb.2012.0021. Epub 2012 Apr 16. PMID: 22506599; PMCID: PMC3342519.
2) Camacho, C. et al 2009. BLAST+: architecture and applications. BMC Bioinformatics, 10, 421.
3) Cock, P.J. et al. 2009. Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics, 25(11), pp.1422–1423.
4) Dinhobl M, Spinard E, Tesler N, Birtley H, Signore A, Ambagala A, Masembe C, Borca MV, Gladue DP. Reclassification of ASFV into 7 Biotypes Using Unsupervised Machine Learning. Viruses. 2023 Dec 30;16(1):67. doi: 10.3390/v16010067. PMID: 38257767; PMCID: PMC10819123.
5) Dinhobl M, Spinard E, Birtley H, Tesler N, Borca MV, Gladue DP. African swine fever virus P72 genotyping tool. Microbiol Resour Announc. 2024 Feb 15;13(2):e0089123. doi: 10.1128/mra.00891-23. Epub 2024 Jan 8. PMID: 38189309; PMCID: PMC10868209.
6) Edgar RC. MUSCLE: multiple sequence alignment with high accuracy and high throughput. Nucleic Acids Res. 2004 Mar 19;32(5):1792-7. doi: 10.1093/nar/gkh340. PMID: 15034147; PMCID: PMC390337.
7) Harris, C.R. et al. Array programming with NumPy. Nature 585, 357–362 (2020). DOI: 10.1038/s41586-020-2649-2.
8) Hunter, J. D. Matplotlib: A 2D Graphics Environment, in Computing in Science & Engineering, vol. 9, no. 3, pp. 90-95, May-June 2007, doi: 10.1109/MCSE.2007.55.
9) Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics, 34:3094-3100. doi:10.1093/bioinformatics/bty191
10) Li, H. (2021). New strategies to improve minimap2 alignment accuracy. Bioinformatics, 37:4572-4574. doi:10.1093/bioinformatics/btab705
11) McKinney, W. et al. 2010. Data structures for statistical computing in python. In Proceedings of the 9th Python in Science Conference. pp. 51–56.
12) Danecek, P. et al. Twelve years of SAMtools and BCFtools, GigaScience, Volume 10, Issue 2, February 2021, giab008, https://doi.org/10.1093/gigascience/giab008
13) Quinlan AR and Hall IM, 2010. BEDTools: a flexible suite of utilities for comparing genomic features. Bioinformatics. 26, 6, pp. 841–842.
14) Chen, S. 2023. Ultrafast one-pass FASTQ data preprocessing, quality control, and deduplication using fastp. iMeta 2: e107. https://doi.org/10.1002/imt2.107
15) Chen, S. et al. ; fastp: an ultra-fast all-in-one FASTQ preprocessor, Bioinformatics, Volume 34, Issue 17, 1 September 2018, Pages i884–i890, https://doi.org/10.1093/bioinformatics/bty560
16) Spinard E, Dinhobl M, Tesler N, Birtley H, Signore AV, Ambagala A, Masembe C, Borca MV, Gladue DP. A Re-Evaluation of African Swine Fever Genotypes Based on p72 Sequences Reveals the Existence of Only Six Distinct p72 Groups. Viruses. 2023 Nov 11;15(11):2246. doi: 10.3390/v15112246. PMID: 38005923; PMCID: PMC10675559.
17) Vasimuddin, Md. et al. Efficient Architecture-Aware Acceleration of BWA-MEM for Multicore Systems. IEEE Parallel and Distributed Processing Symposium (IPDPS), 2019. 10.1109/IPDPS.2019.00041
