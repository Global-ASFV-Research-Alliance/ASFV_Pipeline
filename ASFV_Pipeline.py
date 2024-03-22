import datetime
import sys

from DeNovo_Functions import *

# FILES AND DIRECTORIES
REFERENCEDIRECTORY = "/app/01_References/"                                   #Permanent Directory within Container 
GENOTYPINGREFERENCEFASTAS = "/app/02_AdditionalFiles/p72GenotypingData.csv"  #Permanent File within Container
MEGAGBK = "/app/02_AdditionalFiles/Mega_ASFV_v07.gbk"                        #Permanent File within Container
GRAPHICS = "/app/02_AdditionalFiles/Graphic.txt"                             #Permanent File within Container
BIOTYPINGREFERENCECSV = "/app/02_AdditionalFiles/BiotypingGenomes.csv"       #Permanent File within Container
BIOTYPEREFERENCES = "/app/05_Biotyping_References/"                          #Permanent File within Container
METADATASTORAGE = "MetaDataStorage.csv"                                      #Permanent File Outside of Container, Record of storage
METADATA = "MetadataNew.csv"                                                 #Permanent File Outside of Container, Created from User Input

#For De Novo Pipeline
HOSTREFERENCEDIRECTORY = "/app/06_Sus_scrofa_Reference/"
SUSSCROFAREFERENCE = "Sus_scrofa.Sscrofa11.1.dna.toplevel.fa"

# VARIABLES
ORGANISM = "African swine fever virus"
THREADS = "10"
MinLengthThreshold = 1000
MaxLengthThreshold = 10000
MaxSubSample = 100
Threshold = 0.2
Require_Multiple = 1
Delete_Temp_Files = True

# SCRIPT METADATA                                                                                                                 
TIME = str(datetime.date.today())
VERSION = str("05")
WEBSITE = "Data processed by pipeline (version " + VERSION + ") at " + TIME +  ". This pipeline was constructed by Edward Spinard, Mark Dinhobl, Douglas Gladue, Jacob Fenster, Nicolas Tesler, Hillary Birtley, Cass Erdelyan, James O'Dwyer, and Anthony Signore. Pipeline can be found at https://github.com/Global-ASFV-Research-Alliance/ASFV_Pipeline."

MetaDataAll = MetaDataRead(METADATA)

# START OF CODE..................................................................................................................................................
for name,MetaData in MetaDataAll.iterrows():
    ## Preventing Deletion of Existing Files
    FILESBEFORERUNNING = glob.glob("*")
    
    if MetaData["Project_ID (internal)"] == None:
        ProjectName = TIME + "_" + rand_pass(8) 
        ProjectName = ProjectName.replace(" ", "_").replace("-","_")
    else:
        ProjectName = slugify(MetaData["Project_ID (internal)"]) + "_" + TIME.replace(" ", "_").replace("-","_") + "_" + rand_pass(4)
    MetaData["Project_ID (internal)"] = ProjectName

    MetaDataUpdate(METADATASTORAGE, MetaData)
    print(ProjectName)
    ## Checking Files
    myList = [MetaData["L1R1"], MetaData["L1R2"], MetaData["L2R1"], MetaData["L2R2"], MetaData["L3R1"], MetaData["L3R2"], MetaData["L4R1"], MetaData["L4R2"]]
    myList = res = [i for i in myList if i is not None]
    ErrorString = ""
    
    ## Check that Illumina Directory Exist
    DirectoryExists = True
    for filelocation in myList:
        check_file = os.path.exists(filelocation)
        if check_file == False:
            DirectoryExists = False
    if DirectoryExists == False:
        ErrorString += "Failure. Input Not Found. Make Sure all lanes have a proper name."
    
    ## Check that Minion Directory Exist
    if MetaData["MinionDirectory"] != None:
        if os.path.isdir(MetaData["MinionDirectory"]) == False:
            ErrorString += "Failure. Minion Directory Not Found. Make sure folder is named correctly and with an ending /"
        MetaData["MinionDirectory"] = MetaData["MinionDirectory"] + "/"
    
    ## Check that all lanes are unique
    if len(myList) != len(set(myList)):
        ErrorString += "Failure. Input Repeated. Make Sure all lanes are seperate files."

    if ErrorString != "":
        ErrorUpdate(MetaDataFile=METADATASTORAGE,
            Error=ErrorString,
            ProjectName=ProjectName)
        continue

     
    # Variables of Files to be Created
    BLASTOUTPUTBANK = Make_Directory(ProjectName = ProjectName)
    
    ## Start of Best Match PipeLine.................................................................................................................................
    if MetaData["MinionDirectory"] == None:

        Trim_L1R1, Trim_L1R2, HTML_File1 = TrimReads(Threads=THREADS, ProjectName=ProjectName, Read01 = MetaData["L1R1"], Read02 = MetaData["L1R2"], Lane = 1)
        Trim_L2R1, Trim_L2R2, HTML_File2 = TrimReads(Threads=THREADS, ProjectName=ProjectName, Read01 = MetaData["L2R1"], Read02 = MetaData["L2R2"], Lane = 2)
        Trim_L3R1, Trim_L3R2, HTML_File3 = TrimReads(Threads=THREADS, ProjectName=ProjectName, Read01 = MetaData["L3R1"], Read02 = MetaData["L3R2"], Lane = 3)
        Trim_L4R1, Trim_L4R2, HTML_File4 = TrimReads(Threads=THREADS, ProjectName=ProjectName, Read01 = MetaData["L4R1"], Read02 = MetaData["L4R2"], Lane = 4)

        ## Part 02 Map Reads to Georgia...........................................................................................................................
        MappedToGeorgia_L1 = MapToReference(Trimmed_Read01 = Trim_L1R1,
            Trimmed_Read02 = Trim_L1R2,
            Lane = 1,
            Threads = THREADS,
            ReferenceDirectory = REFERENCEDIRECTORY,
            Reference = "LR743116_Georgia_2007.fa",
            ProjectName = ProjectName)
        MappedToGeorgia_L2 = MapToReference(Trimmed_Read01 = Trim_L2R1,
            Trimmed_Read02 = Trim_L2R2,
            Lane = 2,
            Threads = THREADS,
            ReferenceDirectory = REFERENCEDIRECTORY,
            Reference = "LR743116_Georgia_2007.fa",
            ProjectName = ProjectName)
        MappedToGeorgia_L3 = MapToReference(Trimmed_Read01 = Trim_L3R1,
            Trimmed_Read02 = Trim_L3R2,
            Lane = 3,
            Threads = THREADS,
            ReferenceDirectory = REFERENCEDIRECTORY,
            Reference = "LR743116_Georgia_2007.fa",
            ProjectName = ProjectName)
        MappedToGeorgia_L4 = MapToReference(Trimmed_Read01 = Trim_L4R1,
            Trimmed_Read02 = Trim_L4R2,
            Lane = 4,
            Threads = THREADS,
            ReferenceDirectory = REFERENCEDIRECTORY,
            Reference = "LR743116_Georgia_2007.fa",
            ProjectName = ProjectName)
        # Combine Mappings
        MappedToGeorgia_Combined = MergeMaps(MapL1 = MappedToGeorgia_L1, 
            MapL2 = MappedToGeorgia_L2, 
            MapL3 = MappedToGeorgia_L3, 
            MapL4 = MappedToGeorgia_L4, 
            ReferenceName = "LR743116_Georgia_2007", 
            ProjectName = ProjectName)

        ## Part 03 Extracted Reads Mapped to Georgia Reference & Collect UnMapped Reads.............................................................................
        ExtractedRead01, ExtractedRead02 = ExtractFromMerged(Input = MappedToGeorgia_Combined, Reference = "LR743116_Georgia_2007", ProjectName = ProjectName)
        
        # New 3/14/2024~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        UnMappedGeorgia01, UnMappedGeorgia02 = ExtractFromMerged_UnMappedReads(Input = MappedToGeorgia_Combined, Reference = "UnMapped_Georgia_Extracted", ProjectName = ProjectName)

        ## Part 04 Map Unmapped Georgia Reads to Swine using Stringent Settings & Collect UnMapped Reads............................................................
        Mapped_SusScrofa = MapToReference(Trimmed_Read01 = UnMappedGeorgia01, Trimmed_Read02 = UnMappedGeorgia02, Lane = 1, Threads = THREADS, ReferenceDirectory = HOSTREFERENCEDIRECTORY, Reference = SUSSCROFAREFERENCE, ProjectName = ProjectName, Stringent = True)
        UnMappedSwine01, UnMappedSwine02 = ExtractFromMerged_UnMappedReads(Input = Mapped_SusScrofa, Reference = SUSSCROFAREFERENCE, ProjectName = ProjectName)

        ## Part 05 Spades De Novo...................................................................................................................................
        Assembled_Scaffold_File = SpadesDeNovo(ExtractedRead01, ExtractedRead02, ProjectName, None, UnMappedSwine01, UnMappedSwine02)
        
        # New 3/15/2024~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        # Remove Assembled Swine Contigs and/or rename + move Scaffold file
        Assembled_Scaffold_File = RemoveSwineContigs(DeNovo_Scaffold = Assembled_Scaffold_File, ProjectName = ProjectName)
        
        ## Check if Assembled scaffold file was created
        if Empty_File_Check(Assembled_Scaffold_File) == True:
            ErrorString = "Failure. De Novo Assembly did not create a scaffold file. This is indicative of poor sequencing or a glitch in the SPades software. Please check the following files:" + ExtractedRead01 + ", " + ExtractedRead02  + ", " + UnMappedSwine01  + ", " + UnMappedSwine02 + ". The Pipeline has finished prematurely."
            ErrorUpdate(MetaDataFile=METADATASTORAGE, Error=ErrorString, ProjectName=ProjectName)
            continue

        """
        # Copy and Rename Scaffold file if there are not swine contigs
        if Perform_Rename_Function == True:
            MoveReferenceFile(ReferenceDirectory="" + ProjectName + "_DeNovoAssembly/", PredictedReference = "scaffolds", FileType=".fasta")
            Assembled_Scaffold_File_to_Keep = ProjectName + "_scaffolds.fasta"
            Rename_File(OldFile = "scaffolds.fasta", NewFile = Assembled_Scaffold_File_to_Keep)
        else:
            Assembled_Scaffold_File_to_Keep = Assembled_Scaffold_File
        """
        
        ## Part 06 Predict Reference...................................................................................................................................
        PredictedReference = PredictReference(Assembled_Scaffold_File, ProjectName)
        
        ## Check if PredictedReference Went through
        if PredictedReference == None:
            ErrorString = "Failure. Failed to Predict Reference likely due to a lack of contigs over 1250 in length. Please check your assembly."
            ErrorUpdate(MetaDataFile=METADATASTORAGE, Error=ErrorString, ProjectName=ProjectName)
            continue
        
        ## Part 07a Map Scaffold to Predicted Reference
        ScaffoldMapped_SAM = MiniMapScaffoldToReference(GenomeFile = REFERENCEDIRECTORY + PredictedReference + ".fa",
            ScaffoldFile = Assembled_Scaffold_File, 
            ProjectName = ProjectName, 
            Suffix = "DenovoScaffold_Mapped_to_" + PredictedReference)

        ## Part 07 Map to Predicted Reference
        if PredictedReference == "LR743116_Georgia_2007.fa":
            CombinedPredicted = MappedToGeorgia_Combined
        else:
            PredictedMap_L1 = MapToReference(Trim_L1R1, Trim_L1R2, Lane = 1, Threads = THREADS, ReferenceDirectory = REFERENCEDIRECTORY, Reference=PredictedReference + ".fa", ProjectName=ProjectName)
            PredictedMap_L2 = MapToReference(Trim_L2R1, Trim_L2R2, Lane = 2, Threads = THREADS, ReferenceDirectory = REFERENCEDIRECTORY, Reference=PredictedReference + ".fa", ProjectName=ProjectName)
            PredictedMap_L3 = MapToReference(Trim_L3R1, Trim_L3R2, Lane = 3, Threads = THREADS, ReferenceDirectory = REFERENCEDIRECTORY, Reference=PredictedReference + ".fa", ProjectName=ProjectName)
            PredictedMap_L4 = MapToReference(Trim_L4R1, Trim_L4R2, Lane = 4, Threads = THREADS, ReferenceDirectory = REFERENCEDIRECTORY, Reference=PredictedReference + ".fa", ProjectName=ProjectName)
            CombinedPredicted = MergeMaps(PredictedMap_L1,PredictedMap_L2,PredictedMap_L3,PredictedMap_L4,PredictedReference, ProjectName)


        Ref_File_In_Directory = MoveReferenceFile(ReferenceDirectory=REFERENCEDIRECTORY, PredictedReference=PredictedReference, FileType=".fa")

        print("NewPart")
        Corrected_CallsVCF = Consensus_Extracter(Ref_File_In_Directory, CombinedPredicted, str(1000000), "", ProjectName)

        #Needed!
        UnmappedBed = UnMappedRegions(CombinedPredicted, "", ProjectName)
        ConsensusGenome = Apply_Variants_to_Consensus(Ref_File_In_Directory,Corrected_CallsVCF,UnmappedBed,ProjectName,"consensus")

        ### Right Line Starts Here
        Predicted_Genotype_Info_List = p72finder_new(queryfile = ConsensusGenome, subjectfile = GENOTYPINGREFERENCEFASTAS, ProjectName = ProjectName)
        
        Predicted_Biotype, BioTypeStatement = BiotypingWithBlast(queryfile = ConsensusGenome, 
            blastoutputbank =  BLASTOUTPUTBANK, 
            FastaBank = BIOTYPEREFERENCES, 
            BioTypeReferenceCSV = BIOTYPINGREFERENCECSV)

        ConsensusGenome = MetaDataString(input_fasta_file = ConsensusGenome, 
            isolate = MetaData['isolate'], 
            organism = ORGANISM, 
            collection_date = MetaData['collection_date'], 
            country = MetaData['country'], 
            location = MetaData['location'], 
            host = MetaData['host'], 
            tissue = MetaData['tissue_type'], 
            collected_by = MetaData['collected_by'], 
            Genotype = Predicted_Genotype_Info_List[0],
            Biotype = Predicted_Biotype,
            note = WEBSITE)

        Genbank_File = Annotation_Function(ProjectName = ProjectName, 
            Fasta_to_Annotate = ConsensusGenome, 
            Genbank_to_scrape = MEGAGBK)


        Genome_Index(ConsensusGenome)
        ConsensusMap_L1 = MapToReference(Trim_L1R1, Trim_L1R2, Lane = 1, Threads = THREADS, ReferenceDirectory="", Reference=ConsensusGenome, ProjectName=ProjectName)
        ConsensusMap_L2 = MapToReference(Trim_L2R1, Trim_L2R2, Lane = 2, Threads = THREADS, ReferenceDirectory="", Reference=ConsensusGenome, ProjectName=ProjectName)
        ConsensusMap_L3 = MapToReference(Trim_L3R1, Trim_L3R2, Lane = 3, Threads = THREADS, ReferenceDirectory="", Reference=ConsensusGenome, ProjectName=ProjectName)
        ConsensusMap_L4 = MapToReference(Trim_L4R1, Trim_L4R2, Lane = 4, Threads = THREADS, ReferenceDirectory="", Reference=ConsensusGenome, ProjectName=ProjectName)
        print("ConsensusMap Complete Done")
        MergedMapConsensus = MergeMaps(ConsensusMap_L1,ConsensusMap_L2,ConsensusMap_L3,ConsensusMap_L4,"Consensus", ProjectName=ProjectName)
        Final_Genome_Stats_File = Generate_Mapping_Stats(GenomeFile = MergedMapConsensus, ProjectName= ProjectName)

        ## Part 20 Graph Statistics.................................................................................................................................
        Coverage_Graph_File, Quality_Graph_File, AverageCoverage = Create_Graph(Stats_File = Final_Genome_Stats_File, ProjectName = ProjectName)


        NCBI_Annotation_File = FixGenbankMakeFeatureTable(ProjectName=ProjectName, 
            GBKFile=Genbank_File, 
            ConcensusFasta=ConsensusGenome, 
            statistics_file=Final_Genome_Stats_File)

        Results_File = Create_Results_File(ProjectName = ProjectName, 
            Graphic = GRAPHICS, 
            Website = WEBSITE, 
            isolate = MetaData['isolate'], 
            Time = TIME,
            p72_Genotype = Predicted_Genotype_Info_List[0], 
            p72_Isolate = Predicted_Genotype_Info_List[3], 
            p72_Accession = Predicted_Genotype_Info_List[4], 
            p72_PID = Predicted_Genotype_Info_List[1], 
            p72_Length = Predicted_Genotype_Info_List[2], 
            p72_Warning = Predicted_Genotype_Info_List[5],
            Organism = ORGANISM,
            collection_date = MetaData['collection_date'], 
            country = MetaData['country'], 
            location = MetaData['location'], 
            host = MetaData['host'], 
            tissue = MetaData['tissue_type'], 
            collected_by = MetaData['collected_by'],
            Fasta_Assembly = ConsensusGenome,
            GenBank_Assembly = Genbank_File,
            Biotype = Predicted_Biotype + BioTypeStatement,
            AverageCoverage = AverageCoverage)

        Files_to_Keep_Temp = [HTML_File1, 
            HTML_File2, 
            HTML_File3, 
            HTML_File4, 
            Final_Genome_Stats_File, 
            Coverage_Graph_File, 
            Quality_Graph_File,
            Predicted_Genotype_Info_List[6], 
            ConsensusGenome, 
            Genbank_File, 
            NCBI_Annotation_File,
            MergedMapConsensus, 
            Results_File,
            ScaffoldMapped_SAM,
            Assembled_Scaffold_File,
            CombinedPredicted]
        Files_to_Keep = []
        # Loop to Remove "None" Values
        for file in Files_to_Keep_Temp:
            if file != None :
                Files_to_Keep.append(file)
                
        Move_File(ProjectName = ProjectName, List_of_Files = Files_to_Keep)

    ## Start of DeNovo Pipeline.....................................................................................................................................
    else:
    
        First_Output_Fasta = ProjectName + "_First_10_PartialReads.fa"
        Last_Output_Fasta = ProjectName + "_Last_10_PartialReads.fa"
        First_Output_ALN = ProjectName + "_First_10_Alignment.fa"
        Last_Output_ALN = ProjectName + "_Last_10_Alignment.fa"
        Consensus_with_FlankingFile = ProjectName + "_assembly_with_flanking.fa"

        ## Part 01 Trim Reads .....................................................................................................................................
        Trim_L1R1, Trim_L1R2, HTML_File1 = TrimReads(Threads=THREADS, ProjectName=ProjectName, Read01 = MetaData["L1R1"], Read02 = MetaData["L1R2"], Lane = 1)
        Trim_L2R1, Trim_L2R2, HTML_File2 = TrimReads(Threads=THREADS, ProjectName=ProjectName, Read01 = MetaData["L2R1"], Read02 = MetaData["L2R2"], Lane = 2)
        Trim_L3R1, Trim_L3R2, HTML_File3 = TrimReads(Threads=THREADS, ProjectName=ProjectName, Read01 = MetaData["L3R1"], Read02 = MetaData["L3R2"], Lane = 3)
        Trim_L4R1, Trim_L4R2, HTML_File4 = TrimReads(Threads=THREADS, ProjectName=ProjectName, Read01 = MetaData["L4R1"], Read02 = MetaData["L4R2"], Lane = 4)

        ## Part 02 Map Reads to ASFVG..............................................................................................................................
        MappedToGeorgia_L1 = MapToReference(Trimmed_Read01 = Trim_L1R1,
            Trimmed_Read02 = Trim_L1R2,
            Lane = 1,
            Threads = THREADS,
            ReferenceDirectory = REFERENCEDIRECTORY,
            Reference = "LR743116_Georgia_2007.fa",
            ProjectName = ProjectName)
        MappedToGeorgia_L2 = MapToReference(Trimmed_Read01 = Trim_L2R1,
            Trimmed_Read02 = Trim_L2R2,
            Lane = 2,
            Threads = THREADS,
            ReferenceDirectory = REFERENCEDIRECTORY,
            Reference = "LR743116_Georgia_2007.fa",
            ProjectName = ProjectName)
        MappedToGeorgia_L3 = MapToReference(Trimmed_Read01 = Trim_L3R1,
            Trimmed_Read02 = Trim_L3R2,
            Lane = 3,
            Threads = THREADS,
            ReferenceDirectory = REFERENCEDIRECTORY,
            Reference = "LR743116_Georgia_2007.fa",
            ProjectName = ProjectName)
        MappedToGeorgia_L4 = MapToReference(Trimmed_Read01 = Trim_L4R1,
            Trimmed_Read02 = Trim_L4R2,
            Lane = 4,
            Threads = THREADS,
            ReferenceDirectory = REFERENCEDIRECTORY,
            Reference = "LR743116_Georgia_2007.fa",
            ProjectName = ProjectName)
        # Combine Mappings
        MappedToGeorgia_Combined = MergeMaps(MapL1 = MappedToGeorgia_L1, 
            MapL2 = MappedToGeorgia_L2, 
            MapL3 = MappedToGeorgia_L3, 
            MapL4 = MappedToGeorgia_L4, 
            ReferenceName = "LR743116_Georgia_2007", 
            ProjectName = ProjectName)

        ## Part 03 Extracted Reads Mapped to Georgia Reference & Collect UnMapped Reads.............................................................................
        ExtractedRead01, ExtractedRead02 = ExtractFromMerged(Input = MappedToGeorgia_Combined, Reference = "LR743116_Georgia_2007", ProjectName = ProjectName)
        
        UnMappedGeorgia01, UnMappedGeorgia02 = ExtractFromMerged_UnMappedReads(Input = MappedToGeorgia_Combined, Reference = "UnMapped_Georgia_Extracted", ProjectName = ProjectName)
        
        ## Part 04 Map Unmapped Georgia Reads to Swine using Stringent Settings & Collect UnMapped Reads............................................................
        Mapped_SusScrofa = MapToReference(Trimmed_Read01 = UnMappedGeorgia01, Trimmed_Read02 = UnMappedGeorgia02, Lane = 1, Threads = THREADS, ReferenceDirectory = HOSTREFERENCEDIRECTORY, Reference = SUSSCROFAREFERENCE, ProjectName = ProjectName, Stringent = True)
        UnMappedSwine01, UnMappedSwine02 = ExtractFromMerged_UnMappedReads(Input = Mapped_SusScrofa, Reference = SUSSCROFAREFERENCE, ProjectName = ProjectName)

        ## Part 05 Merge Minion Reads...............................................................................................................................
        # updated 3/13/2024 to get fastq.gz count, needed for file deletion step
        CombinedMinionReadFile, CountofFastqGZ = CombineMinionReads(Directory = MetaData["MinionDirectory"], ProjectName = ProjectName)

        ## Part 06 Spades De Novo...................................................................................................................................
        Assembled_Scaffold_File = SpadesDeNovo(ExtractedRead01 = ExtractedRead01, 
            ExtractedRead02 = ExtractedRead02,
            ProjectName = ProjectName, 
            CombinedMinionReadFile = CombinedMinionReadFile,
            ExtractedRead03 = UnMappedSwine01,
            ExtractedRead04 = UnMappedSwine02
            )

        ## Check if Assembled scaffold file was created
        if Empty_File_Check(Assembled_Scaffold_File) == True:
            ErrorString = "Failure. De Novo Assembly did not create a scaffold file. This is indicative of poor sequencing or a glitch in the SPades software. Please check the following files:" + ExtractedRead01 + ", " + ExtractedRead02  + ", " + UnMappedSwine01  + ", " + UnMappedSwine02 + ". The Pipeline has finished prematurely."
            ErrorUpdate(MetaDataFile=METADATASTORAGE, Error=ErrorString, ProjectName=ProjectName)
            continue

        ## Part 07 Extract Largest Contig from Spades Assembly......................................................................................................
        GenomeFileTemp, Total_Contig_Length = ExtractLargestContig(Assembled_Scaffold_File = Assembled_Scaffold_File, ProjectName = ProjectName)

        # End Script is sequencing quality is poor
        if End_Script(Total_Contig_Length = Total_Contig_Length) == True:
            ErrorString = "Failure. The initial assembly was smaller than 150,000 nt. This is indicative of poor sequencing. The Pipeline has finished prematurely."
            ErrorUpdate(MetaDataFile=METADATASTORAGE,
            Error=ErrorString,
            ProjectName=ProjectName)
            continue
            
        ## Part 08 MiniMap Minion Reads to Largest Contig...........................................................................................................

        # Set Up Variables
        x3_Start_NT = Total_Contig_Length - 10
        x3_Start_NT = str(x3_Start_NT)
        x3_End_NT = str(Total_Contig_Length)

        MiniMap_FileOutputs = MiniMapToReference(GenomeFileTemp = GenomeFileTemp, 
            CombinedMinionReadFile = CombinedMinionReadFile, 
            ProjectName = ProjectName, 
            x5_Start_NT = "0", 
            x5_End_NT = "10", 
            x3_Start_NT = x3_Start_NT, 
            x3_End_NT = x3_End_NT)
        # The List is as follows MiniMap_FileOutputs = [ MinionMappedTempSAM, MinionMappedTempBAM, MinionMappedTempSOR, MinionMappedFirst10, MinionMappedLast10 ]

        ## Part 09 Extract Flanking Sequences, Align Flanking Sequences, Create Consensus, Add Consensus sequences to 5' and 3' end.................................
        x5_Prime_NonPolish_FlankingFile, x3_Prime_NonPolish_FlankingFile = Contig_Extender_for_Polish_V02(ProjectName = ProjectName,
            MinLengthThreshold = MinLengthThreshold,
            MaxLengthThreshold = MaxLengthThreshold,
            MaxSubSample = MaxSubSample, 
            Threshold = Threshold, 
            Require_Multiple = Require_Multiple, 
            MinionMappedFirst10 = MiniMap_FileOutputs[3], 
            MinionMappedLast10 = MiniMap_FileOutputs[4], 
            First_Output_Fasta = First_Output_Fasta, 
            Last_Output_Fasta = Last_Output_Fasta, 
            First_Output_ALN = First_Output_ALN, 
            Last_Output_ALN = Last_Output_ALN)

        ## Part_10 Index the flanking Regions......................................................................................................................
        Genome_Index(Genome_To_Index = x5_Prime_NonPolish_FlankingFile)
        Genome_Index(Genome_To_Index = x3_Prime_NonPolish_FlankingFile)

        ## Part_11 Map Illumina Reads to the flanking regions.......................................................................................................
        Mapped_To_5_Prime_Consensus_L1 = MapToReference(Trimmed_Read01 = Trim_L1R1, 
            Trimmed_Read02=Trim_L1R2, 
            Lane=1, 
            Threads=THREADS, 
            ReferenceDirectory="", 
            Reference=x5_Prime_NonPolish_FlankingFile, 
            ProjectName=ProjectName)
        Mapped_To_5_Prime_Consensus_L2 = MapToReference(Trimmed_Read01 = Trim_L2R1, 
            Trimmed_Read02=Trim_L2R2, 
            Lane=2, 
            Threads=THREADS, 
            ReferenceDirectory="", 
            Reference=x5_Prime_NonPolish_FlankingFile, 
            ProjectName=ProjectName)
        Mapped_To_5_Prime_Consensus_L3 = MapToReference(Trimmed_Read01 = Trim_L3R1, 
            Trimmed_Read02=Trim_L3R2, 
            Lane=3, 
            Threads=THREADS, 
            ReferenceDirectory="", 
            Reference=x5_Prime_NonPolish_FlankingFile, 
            ProjectName=ProjectName)
        Mapped_To_5_Prime_Consensus_L4 = MapToReference(Trimmed_Read01 = Trim_L4R1, 
            Trimmed_Read02=Trim_L4R2, 
            Lane=4, 
            Threads=THREADS, 
            ReferenceDirectory="", 
            Reference=x5_Prime_NonPolish_FlankingFile, 
            ProjectName=ProjectName)

        Mapped_To_3_Prime_Consensus_L1 = MapToReference(Trimmed_Read01 = Trim_L1R1, 
            Trimmed_Read02=Trim_L1R2, 
            Lane=1, 
            Threads=THREADS, 
            ReferenceDirectory="", 
            Reference=x3_Prime_NonPolish_FlankingFile, 
            ProjectName=ProjectName)
        Mapped_To_3_Prime_Consensus_L2 = MapToReference(Trimmed_Read01 = Trim_L2R1, 
            Trimmed_Read02=Trim_L2R2, 
            Lane=2, 
            Threads=THREADS, 
            ReferenceDirectory="", 
            Reference=x3_Prime_NonPolish_FlankingFile, 
            ProjectName=ProjectName)
        Mapped_To_3_Prime_Consensus_L3 = MapToReference(Trimmed_Read01 = Trim_L3R1, 
            Trimmed_Read02=Trim_L3R2, 
            Lane=3, 
            Threads=THREADS, 
            ReferenceDirectory="", 
            Reference=x3_Prime_NonPolish_FlankingFile, 
            ProjectName=ProjectName)
        Mapped_To_3_Prime_Consensus_L4 = MapToReference(Trimmed_Read01 = Trim_L4R1, 
            Trimmed_Read02=Trim_L4R2, 
            Lane=4, 
            Threads=THREADS, 
            ReferenceDirectory="", 
            Reference=x3_Prime_NonPolish_FlankingFile, 
            ProjectName=ProjectName)

        # Combine BAM mapping Files
        Mapped_To_5_Prime_Combined = MergeMaps(MapL1 = Mapped_To_5_Prime_Consensus_L1, 
            MapL2 = Mapped_To_5_Prime_Consensus_L2, 
            MapL3 = Mapped_To_5_Prime_Consensus_L3, 
            MapL4 = Mapped_To_5_Prime_Consensus_L4, 
            ReferenceName = "5_Prime_nonpolished", 
            ProjectName = ProjectName)
        Mapped_To_3_Prime_Combined = MergeMaps(MapL1 = Mapped_To_3_Prime_Consensus_L1, 
            MapL2 = Mapped_To_3_Prime_Consensus_L2, 
            MapL3 = Mapped_To_3_Prime_Consensus_L3, 
            MapL4 = Mapped_To_3_Prime_Consensus_L4, 
            ReferenceName = "3_Prime_nonpolished", 
            ProjectName = ProjectName)

        ## Part_12 Find Variants and Unmapped Regions, apply to flanking region fasta...............................................................................
        x5_Prime_Polished_variant_file = Consensus_Extracter(Genome_Fasta = x5_Prime_NonPolish_FlankingFile, 
            BAM_File = Mapped_To_5_Prime_Combined,
            Max_Depth_Coverage="100",
            Suffix ="x5_prime_flank_consensus",  
            ProjectName = ProjectName)
        x3_Prime_Polished_variant_file = Consensus_Extracter(Genome_Fasta = x3_Prime_NonPolish_FlankingFile, 
            BAM_File = Mapped_To_3_Prime_Combined,
            Max_Depth_Coverage="100",
            Suffix ="x3_prime_flank_consensus", 
            ProjectName = ProjectName)
        #
        UnMappedRegion_5_Prime_file = UnMappedRegions(BAM_File = Mapped_To_5_Prime_Combined, Suffix = "x5_prime_flank", ProjectName = ProjectName)
        UnMappedRegion_3_Prime_file = UnMappedRegions(BAM_File = Mapped_To_3_Prime_Combined, Suffix = "x3_prime_flank", ProjectName = ProjectName)
        #
        x5_Prime_Polished_Flanking = Apply_Variants_to_Consensus(Genome_Fasta = x5_Prime_NonPolish_FlankingFile, 
            VCF_Index = x5_Prime_Polished_variant_file, 
            BED_File = UnMappedRegion_5_Prime_file, 
            ProjectName = ProjectName, 
            Suffix = "x5_prime_flank")
        x3_Prime_Polished_Flanking = Apply_Variants_to_Consensus(Genome_Fasta = x3_Prime_NonPolish_FlankingFile, 
            VCF_Index = x3_Prime_Polished_variant_file, 
            BED_File = UnMappedRegion_3_Prime_file, 
            ProjectName = ProjectName, 
            Suffix = "x3_prime_flank")
        
        # new 3/12/2024
        if Empty_File_Check(x5_Prime_Polished_Flanking) == True:
            x5_Prime_Polished_Flanking = x5_Prime_NonPolish_FlankingFile
        if Empty_File_Check(x3_Prime_Polished_Flanking) == True:
            x3_Prime_Polished_Flanking = x3_Prime_NonPolish_FlankingFile

        ## Part_13 Add flanking Sequence to longest contig..........................................................................................................
        Consensus_with_FlankingFile = Add_Flanking_Sequence(ProjectName = ProjectName, 
            MiddleFasta_File = GenomeFileTemp, 
            x5_PrimeSeqence_File = x5_Prime_Polished_Flanking, 
            x3_PrimeSequence_File = x3_Prime_Polished_Flanking)

        ## Part 14 Index newly created fasta........................................................................................................................
        Genome_Index(Genome_To_Index = Consensus_with_FlankingFile)

        ## Part 15  Map Trimmed Illumina Reads to newly created fasta...............................................................................................
        Mapped_To_Extended_Consensus_L1 = MapToReference(Trimmed_Read01 = Trim_L1R1, 
            Trimmed_Read02=Trim_L1R2, 
            Lane=1, 
            Threads=THREADS, 
            ReferenceDirectory="", 
            Reference=Consensus_with_FlankingFile, 
            ProjectName=ProjectName)
        Mapped_To_Extended_Consensus_L2 = MapToReference(Trimmed_Read01 = Trim_L2R1, 
            Trimmed_Read02=Trim_L2R2, 
            Lane=2, 
            Threads=THREADS, 
            ReferenceDirectory="", 
            Reference=Consensus_with_FlankingFile, 
            ProjectName=ProjectName)
        Mapped_To_Extended_Consensus_L3 = MapToReference(Trimmed_Read01 = Trim_L3R1, 
            Trimmed_Read02=Trim_L3R2, 
            Lane=3, 
            Threads=THREADS, 
            ReferenceDirectory="", 
            Reference=Consensus_with_FlankingFile, 
            ProjectName=ProjectName)
        Mapped_To_Extended_Consensus_L4 = MapToReference(Trimmed_Read01 = Trim_L4R1, 
            Trimmed_Read02=Trim_L4R2, 
            Lane=4, 
            Threads=THREADS, 
            ReferenceDirectory="", 
            Reference=Consensus_with_FlankingFile, 
            ProjectName=ProjectName)
        # Combine BAM mapping Files
        Mapped_To_Extended_Consensus_Combined = MergeMaps(MapL1 = Mapped_To_Extended_Consensus_L1, 
            MapL2 = Mapped_To_Extended_Consensus_L2, 
            MapL3 = Mapped_To_Extended_Consensus_L3, 
            MapL4 = Mapped_To_Extended_Consensus_L4, 
            ReferenceName = "Extended_Consensus", 
            ProjectName = ProjectName)

        ## Part 16 Clean up SNPs present in mapped genome and extract Consensus.....................................................................................
        Final_Genome_preTrim_Variant_File = Consensus_Extracter(Genome_Fasta = Consensus_with_FlankingFile, 
            BAM_File = Mapped_To_Extended_Consensus_Combined,
            Max_Depth_Coverage ="10000",
            Suffix = "PreTrimFinal",
            ProjectName = ProjectName)
        Final_Genome_preTrim_UnMapped_Bed = UnMappedRegions(BAM_File = Mapped_To_Extended_Consensus_Combined, Suffix = "PreTrimFinal", ProjectName = ProjectName)
        Final_Genome_preTrim = Apply_Variants_to_Consensus(Genome_Fasta = Consensus_with_FlankingFile, 
            VCF_Index = Final_Genome_preTrim_Variant_File, 
            BED_File = Final_Genome_preTrim_UnMapped_Bed, 
            ProjectName = ProjectName, 
            Suffix = "PreTrimFinal")

        ## Part 17 Trim Ambgious nucloetides at the end of the genomes..............................................................................................
        Final_Genome_preMeta = Trim_N_at_Start_End(input_fasta = Final_Genome_preTrim, ProjectName = ProjectName)

        ## Part 18 Predicted P72 Genotype...........................................................................................................................
        Predicted_Genotype_Info_List = p72finder_new(queryfile = Final_Genome_preMeta, subjectfile = GENOTYPINGREFERENCEFASTAS, ProjectName = ProjectName)
        """#List = [
        #0 Predicted_Genotype, 
        #1 Predicted_Genotype_PID, 
        #2 Predicted_Genotype_MatchLength, 
        #3 Predicted_Genotype_Match_Isolate, 
        #4 Predicted_Genotype_Match_Accession, 
        #5 WarningString, 
        #6 outputlocation, 
        #7 Predicted_Genotype_qstart, 
        #8 Predicted_Genotype_qend]
        #""" 


        ## Part 19 Reverse Complement Genome if it is in the wrong orientation.....................................................................................
        QStart_Qend_NewFasta(queryfile = Final_Genome_preMeta, qstart = Predicted_Genotype_Info_List[7], qend = Predicted_Genotype_Info_List[8])

        ## Part 20 Index Corrected Fasta............................................................................................................................
        Genome_Index(Genome_To_Index = Final_Genome_preMeta)

        ## Part 21 Map Trimmed Illumina Reads to Corrected Consensus................................................................................................
        Mapped_To_Final_Genome_L1 = MapToReference(Trimmed_Read01 = Trim_L1R1, 
            Trimmed_Read02 = Trim_L1R2, 
            Lane = 1, 
            Threads = THREADS, 
            ReferenceDirectory = "", 
            Reference = Final_Genome_preMeta, 
            ProjectName = ProjectName)
        Mapped_To_Final_Genome_L2 = MapToReference(Trimmed_Read01 = Trim_L2R1, 
            Trimmed_Read02 = Trim_L2R2, 
            Lane = 2, 
            Threads = THREADS, 
            ReferenceDirectory = "", 
            Reference = Final_Genome_preMeta, 
            ProjectName = ProjectName)
        Mapped_To_Final_Genome_L3 = MapToReference(Trimmed_Read01 = Trim_L3R1, 
            Trimmed_Read02 = Trim_L3R2, 
            Lane = 3, 
            Threads = THREADS, 
            ReferenceDirectory = "", 
            Reference = Final_Genome_preMeta, 
            ProjectName = ProjectName)
        Mapped_To_Final_Genome_L4 = MapToReference(Trimmed_Read01 = Trim_L4R1, 
            Trimmed_Read02 = Trim_L4R2, 
            Lane = 4, 
            Threads = THREADS, 
            ReferenceDirectory = "", 
            Reference = Final_Genome_preMeta, 
            ProjectName = ProjectName)

        # Combine BAM mapping Files
        Mapped_To_Final_Genome_Combined = MergeMaps(MapL1 = Mapped_To_Final_Genome_L1, 
            MapL2 = Mapped_To_Final_Genome_L2, 
            MapL3 = Mapped_To_Final_Genome_L3, 
            MapL4 = Mapped_To_Final_Genome_L4, 
            ReferenceName = "Final_Genome", 
            ProjectName = ProjectName)

        ## Part 22 Prep Bam file for analysis, Generate Coverage and Mapping Stats..................................................................................
        Final_Genome_Stats_File = Generate_Mapping_Stats(GenomeFile = Mapped_To_Final_Genome_Combined, ProjectName= ProjectName)

        ## Part 23 Graph Statistics.................................................................................................................................
        Coverage_Graph_File, Quality_Graph_File, AverageCoverage = Create_Graph(Stats_File = Final_Genome_Stats_File, ProjectName = ProjectName)

        ## Part 24 Map Minon Reads to Corrected Consensus (information only)........................................................................................
        MinionMapped_FinalGenome = MiniMapToReferenceSimple2(GenomeFileTemp = Final_Genome_preMeta,
            CombinedMinionReadFile = CombinedMinionReadFile, 
            ProjectName = ProjectName, 
            Suffix = "")

        ## Part 25 Predict Biotype..................................................................................................................................
        Predicted_Biotype, BioTypeStatement = BiotypingWithBlast(queryfile = Final_Genome_preMeta, 
            blastoutputbank =  BLASTOUTPUTBANK, 
            FastaBank = BIOTYPEREFERENCES, 
            BioTypeReferenceCSV = BIOTYPINGREFERENCECSV)

        ## Part 26 Add MetaData to fasta file.........Update with Biotype...........................................................................................
        Final_Genome = MetaDataString(input_fasta_file = Final_Genome_preMeta, 
            isolate = MetaData['isolate'], 
            organism = ORGANISM, 
            collection_date = MetaData['collection_date'], 
            country = MetaData['country'], 
            location = MetaData['location'], 
            host = MetaData['host'], 
            tissue = MetaData['tissue_type'], 
            collected_by = MetaData['collected_by'], 
            Genotype = Predicted_Genotype_Info_List[0],
            Biotype = Predicted_Biotype,
            note = WEBSITE)

        ## Part 27 Annotate Using the Transformer...................................................................................................................
        Genbank_File = Annotation_Function(ProjectName = ProjectName, 
            Fasta_to_Annotate = Final_Genome, 
            Genbank_to_scrape = MEGAGBK)
        """
        ## Part 25 Convert Genbank Annotations to NCBI_annotation file..............................................................................................
        NCBI_Annotation_File = AnnotationToDF(GBKFile = Genbank_File, 
            ProjectName = ProjectName, 
            isolate = MetaData['isolate'])
        """
        NCBI_Annotation_File = FixGenbankMakeFeatureTable(ProjectName=ProjectName, 
            GBKFile=Genbank_File, 
            ConcensusFasta=Final_Genome, 
            statistics_file=Final_Genome_Stats_File)

        ## Part 28 Create Results Text File.........................................................................................................................
        Results_File = Create_Results_File(ProjectName = ProjectName, 
            Graphic = GRAPHICS, 
            Website = WEBSITE, 
            isolate = MetaData['isolate'], 
            Time = TIME,
            p72_Genotype = Predicted_Genotype_Info_List[0], 
            p72_Isolate = Predicted_Genotype_Info_List[3], 
            p72_Accession = Predicted_Genotype_Info_List[4], 
            p72_PID = Predicted_Genotype_Info_List[1], 
            p72_Length = Predicted_Genotype_Info_List[2], 
            p72_Warning = Predicted_Genotype_Info_List[5],
            Organism = ORGANISM,
            collection_date = MetaData['collection_date'], 
            country = MetaData['country'], 
            location = MetaData['location'], 
            host = MetaData['host'], 
            tissue = MetaData['tissue_type'], 
            collected_by = MetaData['collected_by'],
            Fasta_Assembly = Final_Genome,
            GenBank_Assembly = Genbank_File,
            Biotype = Predicted_Biotype + BioTypeStatement,
            AverageCoverage = AverageCoverage)

        #..#..#..#..#..#..#..#..#..#..#..#..#..#..#..#..#..#..#..#..#..#..#..#..#..#..#..#..#..#..#..#..#..#..#..#..#..#..#..#..#..#..#..#..#..#..#..#..#..#..#..#..
        #++ EPILOGUE ++##

        ## Move Files that we want to keep
        Files_to_Keep_Temp = [HTML_File1, 
            HTML_File2, 
            HTML_File3, 
            HTML_File4, 
            Final_Genome_Stats_File, 
            Coverage_Graph_File, 
            Quality_Graph_File,
            Predicted_Genotype_Info_List[6], 
            Final_Genome, 
            Genbank_File, 
            NCBI_Annotation_File, 
            Results_File,
            Mapped_To_Final_Genome_Combined,
            MinionMapped_FinalGenome[2]]
        Files_to_Keep = []
        # Loop to Remove "None" Values
        for file in Files_to_Keep_Temp:
            if file != None :
                Files_to_Keep.append(file)
        Move_File(ProjectName = ProjectName, List_of_Files = Files_to_Keep)

    ## Delete Temp Files............................................................................................................................................
    if Delete_Temp_Files == True:
        files = glob.glob(ProjectName + "*")
        try:
            files += ["fastp.json"]
        except:
            pass
        try:
            files += glob.glob(PredictedReference + "*")
        except:
            pass
        files = [FILE_FOR_DELETION for FILE_FOR_DELETION in files if FILE_FOR_DELETION not in FILESBEFORERUNNING]
        for file in files:
            try:
                Delete_File(file)
            except:
                try: 
                    Delete_Folder(file)
                except: print("Failed to Delete " + file)
        # update 3/13/2024
        try: 
            CountofFastqGZ
            if CountofFastqGZ > 1:
                Delete_File(CombinedMinionReadFile)
        except:
            pass

    ErrorString = ""
    ErrorUpdate(MetaDataFile=METADATASTORAGE,
            Error=ErrorString,
            ProjectName=ProjectName)