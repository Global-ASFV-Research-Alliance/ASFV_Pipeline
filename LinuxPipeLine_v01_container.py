import os
import Bio
import Bio.SeqIO
import pandas as pd
from Bio.Blast.Applications import NcbiblastnCommandline
import time
import glob

import random
import string
import datetime
import matplotlib.pyplot as plt
from statistics import mean
import shutil #NEWLIBRARY

#CSVs
METADATASTORAGE = "MetaDataStorage.csv" #Outside of Container, Record of storage
METADATA = "MetadataNew.csv" #Outside of Container, what they put on the website
GENOTYPINGREFERENCEFASTAS = "p72GenotypingData.csv" #Permanent File
BIOTYPINGREFERENCECSV = "BiotypingGenomes.csv" #Permanent File

#Settings
Space = " "
Threads = "10"
VERSION = 0.2 #Not sure if we are going to do other versions but if we add one that could be inserted here
ReferenceDirectory = "/app/01_References/" #Permanent File
MegaGBK = "/app/TheTransporter-main/Test/Mega_ASFV_v07.gbk" #Permanent File

#New As of 0.2
BiotypingReferences = "/app/Biotyping_References/" #Permanent File

#GivenMetaData
ORGANISM = "African swine fever virus"
WEBSITE = "Data processed by pipeline (version " + str(VERSION) + ") hosted at https://asfvgenomics.com/ on " + str(datetime.date.today()) + "."

#########################################################################

### Setup

#Generate Projectname
def rand_pass(size): 
    # Takes random choices from 
    # ascii_letters and digits 
    generate_pass = ''.join([random.choice( 
                        string.ascii_lowercase + string.digits) 
                        for n in range(size)]) 
                        
    return generate_pass 
ProjectName = str(datetime.date.today()) + "_" + rand_pass(8)
ProjectName = ProjectName.replace(" ", "_").replace("-","_")

#Update MetadataFile
def MetaDataUpdate(MetaDataFile = METADATASTORAGE, NewMetaData = METADATA, ProjectName = ProjectName):
    
    MetaDataFile = pd.read_csv(MetaDataFile)
    MetaDataFile = MetaDataFile.where(pd.notnull(MetaDataFile), None)

    NewMetaData = pd.read_csv(NewMetaData).iloc[0]
    NewMetaData["Project_ID (internal)"] = ProjectName
    
    
    MetaDataFile = pd.concat([MetaDataFile,NewMetaData.to_frame().transpose()],axis=0)
    NewMetaData = NewMetaData.where(pd.notnull(NewMetaData), None)
    
    return MetaDataFile,NewMetaData
MetaDataFile,MetaData = MetaDataUpdate()
MetaDataFile.to_csv(METADATASTORAGE, index=False)

print("MetaDataDone")
#Upload Folder will contain all of an upload. And the metadata. Whenever you upload it you will an identifier to each file (not the metadata) to indicate lane and direction
#For example you might rename a file from marysue.fastq.gz to marysue_lane01_direction_01.fastq.gz

#########################################################################
#1st Loop
#ReadDirectory = "/Users/jacobfenster/Local_Documents/Linux_container_fastq/" #this is the only path specific to the users PC #I'm not sure if we still need this - the server should produce a code.

#### TrimReads Function
def TrimReads(Read01, Read02, Lane = 1, Trim_5_Prim = "25", Trim_3_Prim = "5", Ambig = "2", MinLength = "50", QualityPhred = "20", TrimmedRead01Output = None, TrimmedRead02Output = None):
    
    if Read01 == None or Read02 == None:
        Trimmed_Read01 = None
        Trimmed_Read02 = None
        return Trimmed_Read01, Trimmed_Read02
    
    #Setting Output Names
    if TrimmedRead01Output == None:
        Trimmed_Read01 = ProjectName + "_Lane" + str(Lane) + "_Trimmed1.fastq.gz"
        #Trimmed_Read01 = ProjectName + Read01.replace(".fastq","_Trimmed.fastq")
    else:
        Trimmed_Read01 = TrimmedRead01Output
    if TrimmedRead02Output == None:
        Trimmed_Read02 = ProjectName + "_Lane" + str(Lane) + "_Trimmed2.fastq.gz"
        #Trimmed_Read02 = ProjectName + Read02.replace(".fastq","_Trimmed.fastq")
    else:
        Trimmed_Read02 = TrimmedRead02Output
    
    os.system ("fastp -V --detect_adapter_for_pe -i " + Read01 + " -I " + Read02 + " -o " + Trimmed_Read01 + " -O " + Trimmed_Read02 + " -f " + Trim_5_Prim + " -F " + Trim_5_Prim + " -t " + Trim_3_Prim + " -T " + Trim_3_Prim + " -n " + Ambig + " -l " + MinLength + " -w " + Threads + " -q " + QualityPhred + " --html " + ProjectName + "_Lane" + str(Lane) + ".html")
    return Trimmed_Read01, Trimmed_Read02

Trim_L1R1, Trim_L1R2 = TrimReads(MetaData["L1R1"], MetaData["L1R2"], Lane = 1)
Trim_L2R1, Trim_L2R2 = TrimReads(MetaData["L2R1"], MetaData["L2R2"], Lane = 2)
Trim_L3R1, Trim_L3R2 = TrimReads(MetaData["L3R1"], MetaData["L3R2"], Lane = 3)
Trim_L4R1, Trim_L4R2 = TrimReads(MetaData["L4R1"], MetaData["L4R2"], Lane = 4)
print("TrimDone")
### Map to Reference Function
def MapToReference(Trimmed_Read01, Trimmed_Read02, Lane = 1, Threads = Threads, ReferenceDirectory = ReferenceDirectory, Reference = "LR743116_Georgia_2007.fa", ProjectName = ProjectName):
    #docker container directory
    # JF - this is giving error "ERROR! Unable to open the file: /app/01_References/Georgia-2007-LR743116.fa.bwt.2bit.64"
    #RunCode and Pipe Into SamTools to compress sam to bam
    
    if Trimmed_Read01 == None or Trimmed_Read02 == None:
        OutputLocation = ""
        return OutputLocation
    
    ReferenceLocation = ReferenceDirectory + Reference
    OutputLocation = ProjectName + "_mapped_" + Reference.replace(".fa","") + "_Lane" + str(Lane) + ".bam"
    os.system("bwa-mem2 mem -t " + Threads + " " + ReferenceLocation + " " + Trimmed_Read01 + " " + Trimmed_Read02 + " | samtools sort -o " + OutputLocation + " " + "-")
    return OutputLocation

GeorgiaMap_L1 = MapToReference(Trim_L1R1, Trim_L1R2, Lane = 1)
GeorgiaMap_L2 = MapToReference(Trim_L2R1, Trim_L2R2, Lane = 2)
GeorgiaMap_L3 = MapToReference(Trim_L3R1, Trim_L3R2, Lane = 3)
GeorgiaMap_L4 = MapToReference(Trim_L4R1, Trim_L4R2, Lane = 4)
#time.sleep(60)
print("MapToGeorgiaDone")

#Combine Lanes, or if only lane 1 set only to that lane.

def MergeMaps(MapL1, MapL2, MapL3, MapL4, ReferenceName, ProjectName = ProjectName):
    if str(MapL2) + str(MapL3) + str(MapL4) == "":
        CombinedMap = MapL1
    else:
        CombinedMap = ProjectName + "_MergedMap_" + ReferenceName + ".bam"
        Command = "samtools merge " + CombinedMap + " " + str(MapL1) + " " +  str(MapL2) + " " + str(MapL3) + " " + str(MapL4)
        os.system(Command.replace("  ", " ").replace("  ", " ").replace("  ", " "))
    return CombinedMap

CombinedGeorgia = MergeMaps(GeorgiaMap_L1, GeorgiaMap_L2, GeorgiaMap_L3, GeorgiaMap_L4, "LR743116_Georgia_2007")
print("MergedGeorgia Done")

def ExtractFromMerged(Input, Reference, ProjectName = ProjectName):
    
    os.system("samtools view -b -f 0x2 " + Input + " | samtools sort -n -o " + ProjectName + "_mappedReads_" + Reference + "_sort.bam")
    
    #os.system("samtools view -b -f 0x2 " + Input + " > " + ProjectName + "_mappedReads_" + Reference + "_ext.bam")
    #os.system("samtools sort -n " + ProjectName + "_mappedReads_" + Reference + "_ext.bam -o " + ProjectName + "_mappedReads_" + Reference + "_sort.bam")
    
    #extract as fastq-paired Reads
    ExtractedRead01 = ProjectName + "_" + Reference + "_TrimmedExtracted1.fastq.gz"
    ExtractedRead02 = ProjectName + "_" + Reference + "_TrimmedExtracted2.fastq.gz"
    os.system("samtools fastq -@ 8 " + ProjectName + "_mappedReads_" + Reference + "_sort.bam -1 " + ExtractedRead01 + " -2 " + ExtractedRead02 + " -0 /dev/null -s /dev/null -n")
    return ExtractedRead01, ExtractedRead02

ExtractedGeorgia01, ExtractedGeorgia02 = ExtractFromMerged(Input = CombinedGeorgia, Reference = "LR743116_Georgia_2007")

print("Extracted Done")

def SpadesDeNovo(ExtractedRead01, ExtractedRead02, ProjectName = ProjectName, threads = "3", RAM = "10", Kmers = "21,33,55,75,99"):
    #RAM is in Gigabytes
    OutputDirectory = ProjectName + "_DeNovoAssembly"
    os.system("spades.py --careful -1 " + ExtractedRead01 + " -2 " + ExtractedRead02 + " -t " + str(threads) + " -m " + str(RAM) + " -k " + str(Kmers) + " -o " + OutputDirectory)
    return OutputDirectory

DeNovoGeorgia = SpadesDeNovo(ExtractedGeorgia01, ExtractedGeorgia02)

print("Spades Done")

def PredictReference(DeNovo, ProjectName = ProjectName):
    BlastInput =  DeNovo + "/scaffolds.fasta"
    BlastReference = "/app/01_References/ReferenceFastas.fa" #docker container file
    BlastOutput = ProjectName + "_predictreference_blastoutput.csv"
    os.system(str(NcbiblastnCommandline(cmd='blastn', query = BlastInput, subject= BlastReference, max_hsps = 1, out = BlastOutput, outfmt = "6 qseqid sseqid evalue pident bitscore length qlen slen qstart qend sstart send", evalue=0.001)))
    PredictedReference = pd.read_csv(BlastOutput, sep = '\t', names=('qseqid', 'sseqid', 'evalue', 'pident', 'bitscore', 'length', 'qlen', 'slen', 'qstart', 'qend', 'sstart', 'send')).groupby('sseqid')['bitscore'].sum().sort_values(ascending=False).reset_index().iloc[0]['sseqid'] 
    return PredictedReference

PredictedReference = PredictReference(DeNovo=DeNovoGeorgia) 

print("PredictReference Done")

#Map to Predicted References
if PredictedReference == "LR743116_Georgia_2007.fa":
    #PredictedMap_L1 = GeorgiaMap_L1
    #PredictedMap_L2 = GeorgiaMap_L2
    #PredictedMap_L3 = GeorgiaMap_L3
    #PredictedMap_L4 = GeorgiaMap_L4
    CombinedPredicted = CombinedGeorgia
    #ExtractedPredicted01 = ExtractedGeorgia01
    #ExtractedPredicted02 = ExtractedGeorgia02
    #DeNovoPredicted = DeNovoGeorgia #Do we repeat this?
else:
    PredictedMap_L1 = MapToReference(Trim_L1R1, Trim_L1R2, Lane = 1, Reference=PredictedReference + ".fa")
    PredictedMap_L2 = MapToReference(Trim_L2R1, Trim_L2R2, Lane = 2, Reference=PredictedReference + ".fa")
    PredictedMap_L3 = MapToReference(Trim_L3R1, Trim_L3R2, Lane = 3, Reference=PredictedReference + ".fa")
    PredictedMap_L4 = MapToReference(Trim_L4R1, Trim_L4R2, Lane = 4, Reference=PredictedReference + ".fa")
    #time.sleep(60)
    CombinedPredicted = MergeMaps(PredictedMap_L1,PredictedMap_L2,PredictedMap_L3,PredictedMap_L4,PredictedReference)
    

print(CombinedPredicted)

print("MergePredict Done")


def MoveReferenceFile(ReferenceDirectory, PredictedReference, FileType):
    if not os.path.isfile(PredictedReference + FileType):
        shutil.copy2(ReferenceDirectory + PredictedReference + FileType, os.getcwd())
    else:
        pass
    return PredictedReference + FileType
#

Ref_File_In_Directory = MoveReferenceFile(ReferenceDirectory=ReferenceDirectory, PredictedReference=PredictedReference, FileType=".fa")
### Extract Consensus Sequence and Create Fasta
#old....os.system("bcftools mpileup -Ou -f "+ ReferenceDirectory + PredictedReference + ".fa" + " " + CombinedPredicted + " -d 1000000 | bcftools call -mv -Ov -o "  + ProjectName + "_calls.vcf")
os.system("bcftools mpileup -Ou -f " + Ref_File_In_Directory + " " + CombinedPredicted + " -d 1000000 | bcftools call -mv -Ov -o "  + ProjectName + "_calls.vcf")

RAW_File = ProjectName + "_calls.vcf"
CSV_Output = ProjectName + "_SNP.csv"
VCF_Output = ProjectName + "_corrected_calls.vcf"

###############################################################################

#the first 27 rows can be skipped because they are not formatted properly
Variant_Table = pd.read_csv(RAW_File, sep="\t", skiprows = 27)
RAW_Table = pd.read_csv(RAW_File, sep="\t", skiprows = 27)
#the INFO column has an incredible amount of data that is seperated by ';'
Variant_Table['INFO'] = Variant_Table['INFO'].str.split(';').fillna(Variant_Table['INFO'])
#explode the data into multiple rows but keep the original index # on each created row
Variant_Table=Variant_Table.explode('INFO',ignore_index=False)
#of the data in info, we only care about 'DP' (total coverage) and DP4 (Ref Read F, Ref Read R, Variant Read F, Variant Read R)
searchfor = ['DP']
Variant_Table = Variant_Table[Variant_Table.INFO.str.contains('|'.join(searchfor))]
#Create another dataframe where DP and DP4 is seperated using '=', name the columns ID and Value
SplitColumn = Variant_Table['INFO'].str.split(pat = '=', expand = True)
SplitColumn.columns = ['ID', 'Value']
#Create DP and DP4 columns, since indeces were preserved, the values for DP and DP4 at the same position will be preserved
SplitColumn = SplitColumn.pivot(columns = 'ID', values = 'Value')
#Restore string to integer
CoverageColumn = SplitColumn['DP'].astype(int)
#since DP4 is composed of (Ref Read F, Ref Read R, Variant Read F, Variant Read R), we only care about Variant Read F & R
DP4Column = SplitColumn['DP4'].str.split(',', expand = True).astype(int)
#The sum of this value gives us the number of reads that contain the variant
VariantSum = DP4Column[2] + DP4Column[3]
#The frequency of reads that compare the variants as compared to the total coverage at that position
Frequency = round(VariantSum / CoverageColumn, 3)
 
#Clean Up Table
DropColumns = ['ID', 'FILTER', 'INFO', 'FORMAT', CombinedPredicted]
Variant_Table = Variant_Table[~Variant_Table.index.duplicated(keep='first')].drop(DropColumns, axis = 1)
Variant_Table.rename(columns={
    "#CHROM": "Reference_Strain",
    "POS": "Position",
    "REF": "Reference",
    "ALT": "Sequenced_Sample"
    }
)
#Create columns with calculated values
Variant_Table['Count'] = VariantSum
Variant_Table['Coverage'] = CoverageColumn
Variant_Table['Frequency'] = Frequency
 
#Drop Variants less than 50%
Variant_Table = Variant_Table[Variant_Table.Frequency >= 0.50]
 
#Drop indices in orginal file that are not in the final SNP Table
RAW_Table = RAW_Table[RAW_Table.index.isin(Variant_Table.index)]

#Save Variant File in an User Friendly Format
Variant_Table.to_csv(CSV_Output, sep=',', index=False)
 
#Creating A Corrected Variant File for PipeLine
with open(RAW_File) as myfile:
    first_27_lines = myfile.readlines()[0:27]
myfile.close()
with open(VCF_Output,'a') as Corrected_File:
    for line in first_27_lines:
        Corrected_File.write(line)
Corrected_File.close()
RAW_Table.to_csv(VCF_Output, sep = '\t', mode='a', index = False, header=True)



print("SNP Table Created")

# Create a BedFile Marking Unmapped Regions
os.system("bedtools genomecov -ibam " + CombinedPredicted + " -bga | grep -w 0$ > " + ProjectName + "_unmappedRegions.bed")

print("bedtools genomecov Done")
os.system("bcftools view "+ VCF_Output + " -Oz --write-index -o " + VCF_Output + ".gz")
# apply variants to create consensus sequence
ConsensusGenome = ProjectName + "_consensus.fa"
os.system("cat " + ReferenceDirectory + PredictedReference + ".fa | bcftools consensus " + VCF_Output + ".gz" + " -m " + ProjectName + "_unmappedRegions.bed > " + ConsensusGenome)
print("Consensus Genome Done")  

######################################################################### RIGHT LINE STARTS HERE #########################################################################
outputname_fasta = ProjectName + '_P72.fasta'
blastxout = ProjectName + '_blastxout.txt'
def p72finder(queryfile, subjectfile, GenomeColumn = 'Genome', SequenceColumn = 'Sequence', IsolateColumn = 'Isolate', GenotypeColumn = 'GenotypeNumber', outputname_fasta = outputname_fasta, outputlocation = blastxout):
    
    subject_excel = pd.read_csv(subjectfile, usecols=[GenomeColumn,SequenceColumn,IsolateColumn,GenotypeColumn])
    #subject_excel = pd.read_excel(subjectfile, usecols=[GenomeColumn,SequenceColumn,IsolateColumn,GenotypeColumn])
    with open(outputname_fasta, "w") as f:
        pass

    for index, row in subject_excel.iterrows():      
        with open(outputname_fasta, "a") as f:
            print(">" + row[GenomeColumn] + "\n" + str(row[SequenceColumn]).replace("[", "").replace("]", "").replace("'", ""), file = f)

    #Identify if Nucleotide or Ammino Acid Sequence
    
    for record in Bio.SeqIO.parse(queryfile, "fasta"):
        ok = r'ACTGUWSMKRYBDHVN*'
        Type = all(c in ok for c in str(record.seq))
    
    import os
    if Type == True:
        from Bio.Blast.Applications import NcbiblastxCommandline
        blastcmdline = NcbiblastxCommandline(cmd='blastx', query = queryfile, subject= outputname_fasta, max_hsps = 1, max_target_seqs = 43, out = outputlocation, outfmt = "6 qseqid sseqid evalue pident bitscore length qlen slen qstart qend sstart send qseq sseq", evalue=0.001)
        os.system(str(blastcmdline))
    if Type == False:
        from Bio.Blast.Applications import NcbiblastpCommandline
        blastcmdline = NcbiblastpCommandline(cmd='blastp', query = queryfile, subject= outputname_fasta, max_hsps = 1, max_target_seqs = 43, out = outputlocation, outfmt = "6 qseqid sseqid evalue pident bitscore length qlen slen qstart qend sstart send qseq sseq", evalue=0.001)
        os.system(str(blastcmdline))

    #ReadBlastResults
    prediction = pd.read_csv(outputlocation, sep = '\t', names=('qseqid', 'sseqid', 'evalue', 'pident', 'bitscore', 'length', 'qlen', 'slen', 'qstart', 'qend', 'sstart', 'send', 'qseq', 'sseq')).drop_duplicates(["qseq"])
    result = prediction.merge(subject_excel, left_on='sseqid', right_on=GenomeColumn)
    
    MaxBitscore = result['bitscore'].max()
    MaxPident = result[result['bitscore'] == MaxBitscore]['pident'].max()

    # print(MaxPident)
    result = result[(result['bitscore'] == MaxBitscore) & (result['pident'] == MaxPident)]
        
    
    print('Predicted Genotype: ' + str(result.iloc[0][GenotypeColumn]) + "\n" + "The B646L(P72) encoded by your genome is " + str(result.iloc[0]['pident']) + "%" + " identical to the " + str(result.iloc[0]['length']) + ' ammino acids in ' + str(result.iloc[0][IsolateColumn]) + ' (' + str(result.iloc[0][GenomeColumn]) + ')')
    
    if result.iloc[0]['length'] <= 600:
        print("WARNING!!! No full length B646L found in the submitted sequence. A manual double-check of raw blast output is recommended, and genotyping results are highly likely to be illegitamite.")
    elif result.iloc[0]['length'] <= 635:
        print("Warning - No full length B646L found in submitted sequence. A manual double-check of the raw blast output is recommended, there are potentially multiple genotypes overlapping.")

    return str(result.iloc[0][GenotypeColumn])

Genotype = p72finder(queryfile = ConsensusGenome, subjectfile = GENOTYPINGREFERENCEFASTAS) #Tried "/app/" + GENOTYPINGREFERENCEFASTAS, "/app/" + ConsensusGenome
print("Genotype Done")


def FolderPathFixer(FolderPath):
    """
    An internal function used to fix folder paths.
    
    Parameters
    ---
    A path to a folder.
    
    Returns
    ---
    The folder path, if not made, will be made. Additionally, a "\\" will be added at the end to allow more items to the end.
    """
    import os
    os.makedirs(FolderPath, exist_ok=True)
    FolderPath = FolderPath + "\\"
    return FolderPath.replace("\\\\","\\")

## Create a new directory
def Make_Directory(ProjectName, Suffix):
    NewDirectory = ProjectName + "_" + Suffix +"/"
    if not os.path.exists(NewDirectory):
        os.makedirs(NewDirectory)
        return NewDirectory
    else:
        return NewDirectory

blastoutputbank = Make_Directory(ProjectName = ProjectName, Suffix="BLASTOUTPUT")
fastabank = FolderPathFixer(BiotypingReferences)
print(fastabank)

def BiotypingWithBlast(queryfile, blastoutputbank, FastaBank, GenomeCsv, SaveTempFiles = False):
    
    print(queryfile)
    
    for record in Bio.SeqIO.parse(queryfile, "fasta"):
        ok = r'ACTGUWSMKRYBDHVN*'
        Type = all(c in ok for c in str(record.seq))
    import glob
    import os
    my_files = glob.glob(blastoutputbank + '*')
    for file in my_files:
        os.remove(file)

    my_files = glob.glob(FastaBank + '*')
    for file in my_files:
        filename = file.replace(FastaBank, '').replace('.fasta', '')
        outputlocation = blastoutputbank + filename + '.txt'
        #print(filename)
        import os
        if Type == True:
            from Bio.Blast.Applications import NcbiblastxCommandline
            blastcmdline = NcbiblastxCommandline(cmd='blastx', query = queryfile, subject= file, max_hsps = 1, out = outputlocation, outfmt = "6 qseqid sseqid evalue pident bitscore length qlen slen qstart qend sstart send qseq sseq", evalue=0.001)
            os.system(str(blastcmdline))
        if Type == False:
            from Bio.Blast.Applications import NcbiblastpCommandline
            blastcmdline = NcbiblastpCommandline(cmd='blastp', query = queryfile, subject= file, max_hsps = 1, out = outputlocation, outfmt = "6 qseqid sseqid evalue pident bitscore length qlen slen qstart qend sstart send qseq sseq", evalue=0.001)
            os.system(str(blastcmdline))

    def BlastBankToPidentMatrix(blastoutputbank):
        import glob
        import pandas as pd
        my_files = glob.glob(blastoutputbank + '*')
        DF = pd.DataFrame()
        DF2 = pd.DataFrame()
        DF3 = pd.DataFrame()
        for file in my_files:
            filename = file.replace(blastoutputbank, '').replace('.txt', '')
            prediction = pd.read_csv(file, sep = '\t', names=('qseqid', 'sseqid', 'evalue', 'pident', 'bitscore', 'length', 'qlen', 'slen', 'qstart', 'qend', 'sstart', 'send', 'qseq', 'sseq'))
            for name, pred in prediction.iterrows():
                DF.loc[pred['sseqid'],[filename]] = pred['pident']/100
                if pred['length']/pred['slen'] > 1:
                    Ratio = 1
                else:
                    Ratio = pred['length']/pred['slen']
                DF2.loc[pred['sseqid'],[filename]] = Ratio
                DF3.loc[pred['sseqid'],[filename]] = ((pred['pident']/100)**2)*Ratio
        return DF, DF2, DF3

    BM_Pident, BM_Ratio, BM_PidentRatio = BlastBankToPidentMatrix(blastoutputbank)
    
    if SaveTempFiles == True:
        BM_Pident.to_csv(ProjectName + "_TempPident.csv")
        BM_Ratio.to_csv(ProjectName + "_TempRatio.csv")
        BM_PidentRatio.to_csv(ProjectName + "_TempPidentRatio.csv")
    
    BM_Pident_Mean = BM_Pident.mean(axis=1, numeric_only = True)
    BM_Ratio_Mean = BM_Ratio.mean(axis=1, numeric_only = True)
    BM_PidentRatio_Mean = BM_PidentRatio.mean(axis=1, numeric_only = True)
    BM_Count = BM_PidentRatio.count(axis=1, numeric_only= True)

    Score = (BM_PidentRatio_Mean*(BM_Count**0.5)).sort_values(ascending=False)
    
    if SaveTempFiles == True:
        Score.to_csv(ProjectName + "_TempScore.csv")
    
    B646MatchPercentage = BM_Pident['B646L'][GenomeCsv.loc[Score.index[0]].name]*100

    if B646MatchPercentage >= 0.950:
        if B646MatchPercentage < 0.980:
            print("Warning! B646 Match Percentage is beneath 0.98. This could mean there is potential poor sequencing")
            + "\n"
        print(
            "This software finds the closest matching Biotype through a gene-by-gene match, meaning it matches the data submitted against all ASFV genes in our database."
            + "\n"
            + "Predicted Biotype: " + str(GenomeCsv.loc[Score.index[0]]['Biotype'])
            + "\n"
            + "Weighted Percent Gene-By-Gene Match (Closest in Biotype): " + str(BM_PidentRatio_Mean[Score.index[0]]*100) 
            + "\n"
            + "Percent B646L (P72) Match (Closest in Biotype): " + str(B646MatchPercentage)
            + "\n"
            + "Genes Found (Closest in Biotype): " + str(BM_Count[Score.index[0]]) 
            )
        if BM_PidentRatio_Mean[Score.index[0]] <= 0.975:
            Statement = " Based on weighted percent identity, the genome exhibited a" + str(BM_PidentRatio_Mean[Score.index[0]]*100) + "% similarity to its nearest match within the specified Biotype. This percent identity is below our recommended cutoff of 97.5%. Nevertheless, this outcome does not imply that the genome signifies a novel Biotype, however this genome should be analyzed using the complete algorithm described in 'Reclassification of ASFV into 7 Biotypes Using Unsupervised Machine Learning' (Dinhobl et al 2023). Please contact administration for assistance."
        else:
            Statement = ""
        return str(GenomeCsv.loc[Score.index[0]]['Biotype']) + Statement
    else: 
        return "Warning! B646 match percentage is " + str(B646MatchPercentage) + " which is very low. Ensure proper sequencing and cleaning before submitting."      

Biotype = BiotypingWithBlast(queryfile = ConsensusGenome, blastoutputbank =  blastoutputbank, FastaBank = BiotypingReferences, GenomeCsv = pd.read_csv(BIOTYPINGREFERENCECSV, index_col='Accession'))

print("Biotyping Complete!")

def MetaDataString(organism, collection_date, country, location, host, tissue, collected_by, isolate, Genotype, Biotype, note):
    string = str()
    if isolate != None:
        string += isolate + " [Isolate=" + isolate +"]"
    if organism != None:
        string += " [Organism=" + organism + "]"
    if collection_date != None:
        string += " [Collection_date=" + str(collection_date) + "]"
    if country != None:
        if location != None:
            string += " [Country=" + country + ": " + location + "]"
        else:
            string += " [Country=" + country + "]"
    if host != None:
        string += " [Host=" + host + "]"
    if tissue != None:
        string += " [Tissue_type=" + tissue + "]"
    if collected_by != None:
        string += " [Collected_by=" + collected_by + "]"
    if str(Genotype) != None:
        string += " [Genotype=" + str(Genotype) + "]"
    if str(Biotype) != None:
        string += " [Biotype=" + str(Biotype) + "]"
    if note != None:
        string += " [Note=" + note+"]"
    return string

Note = WEBSITE #In Case We Want to Add More to It Later
MetaDataAttributes = MetaDataString(ORGANISM, MetaData['collection_date'],MetaData['country'],MetaData['location'],MetaData['host'],MetaData['tissue_type'],MetaData['collected_by'], MetaData['isolate'], Genotype, Biotype, Note)

FastaFile = list(Bio.SeqIO.parse(ConsensusGenome, format = 'fasta'))
ConsensusGenomeMetaData = ConsensusGenome.replace('.fa', '_metadata.fa')
with open(ConsensusGenomeMetaData, "w") as f:
    print(">" + MetaData['isolate'] + MetaDataAttributes, file = f)
    print(str(FastaFile[0].seq), file = f) 
print("AddMetaData Complete")

#Annotation
AnnotatedGenbank = ProjectName + "_Annotated"
os.system("TheTransporter -t" + ConsensusGenome + " -a " + MegaGBK + " -o " + AnnotatedGenbank)
print("TheTransporterComplete")

#AnnotationToDF
def AnnotationToDF(GBKFile, Name):
    GBKFile = Bio.SeqIO.read(open(GBKFile, "r"), "genbank")

    dfFeatures = pd.DataFrame(columns=['Isolate', 'Gene','Sequence', 'Start', 'End', 'Strand']) #Isolate Gene Sequence Start End
    for Feature in GBKFile.features:
        
        #Only Check CDS
        if Feature.type != 'CDS':
            continue
        #Gene Name
        try:
            FeatureName = Feature.qualifiers["label"][0]
        except:
            try:
                FeatureName = Feature.qualifiers["gene"][0]
            except:
                try: 
                    FeatureName = Feature.qualifiers["name"][0]
                except:
                    FeatureName = Feature.qualifiers["product"][0]
        
        FeatureSequence = Feature.qualifiers["translation"][0]
        FeatureStart = int(Feature.location.start) + 1
        FeatureEnd = int(Feature.location.end)
        FeatureStrand = Feature.strand
        
        #Construct the DataFrame
        dfFeatures.loc[len(dfFeatures.index)] = [Name,FeatureName,FeatureSequence,FeatureStart,FeatureEnd,FeatureStrand]
    dfFeatures = dfFeatures[dfFeatures['Gene'] != 'putative']
    return dfFeatures

#AnnotationCSV
annotationdf = AnnotationToDF(AnnotatedGenbank + ".gb", ProjectName)
annotationdf.to_csv(ProjectName+"_AnnotationTable.csv", index=False)

with open(ProjectName + "_FeatureTable.txt", "w") as f:
    print(">Feature " + MetaData['isolate'], file = f)
    for name, row in annotationdf.iterrows():
        if row['Strand'] == 1:
            print(str(row['Start']) + "\t" + str(row['End']) + "\t" + "CDS", file=f)
        if row['Strand'] == -1:
            print(str(row['End']) + "\t" + str(row['Start']) + "\t" + "CDS", file=f)
        print("\t" + "product " + str(row['Gene']), file=f)

######################################################################### BOTTOM LINE STARTS HERE #########################################################################

os.system("bwa-mem2 index " + ConsensusGenome)  
print("IndexConsensus Complete")

ConsensusMap_L1 = MapToReference(Trim_L1R1, Trim_L1R2, Lane = 1, ReferenceDirectory="", Reference=ConsensusGenome)
ConsensusMap_L2 = MapToReference(Trim_L2R1, Trim_L2R2, Lane = 2, ReferenceDirectory="", Reference=ConsensusGenome)
ConsensusMap_L3 = MapToReference(Trim_L3R1, Trim_L3R2, Lane = 3, ReferenceDirectory="", Reference=ConsensusGenome)
ConsensusMap_L4 = MapToReference(Trim_L4R1, Trim_L4R2, Lane = 4, ReferenceDirectory="", Reference=ConsensusGenome)
#time.sleep(60)
print("ConsensusMap Complete Done")
MergedMapConsensus = MergeMaps(ConsensusMap_L1,ConsensusMap_L2,ConsensusMap_L3,ConsensusMap_L4,"Consensus")
print("MergedMapConsensus Done")
os.system("samtools view -b -f 0x2 " + MergedMapConsensus +" | samtools sort -o " + ProjectName + "_merge_map_consensus_sort.bam") #-n flag will cause mpileup to fail
print("Consensus View")
os.system("samtools mpileup " + ProjectName + "_merge_map_consensus_sort.bam" + " -s -a | cut -f3,5,6 --complement > " + ProjectName + "_statistics.txt")

print("Consensus samtools mpileup")

###Graph Statistics

plt.close("all")
StatDataFrame = pd.read_csv(ProjectName + "_statistics.txt", sep='\t', header = None, names = ['Genome', 'Position', 'Coverage', 'Quality'], quoting=3)
#Graph Coverage and Save PNG
CoverageGraph = StatDataFrame.plot(kind = 'line',
    x = 'Position',
    y = 'Coverage',
    logy = True,
    legend = False,
    title = 'Read Coverage',
    ylabel = 'Depth of Coverage',
    figsize=(20, 10)
)
CoverageFigure = CoverageGraph.get_figure()
CoverageFigure.savefig(ProjectName + '_CoverageGraph.png')
plt.close("all")
print("CoverageFigure Done")

#Graph Mapping Quality and Save PNG
StatDataFrame['Quality'] = StatDataFrame['Quality'].apply(lambda x: round(mean([ord(i) for i in x]), 1) - 33)
QualityGraph = StatDataFrame.plot(kind = 'line',
    x = 'Position',
    y = 'Quality',
    legend = False,
    title = 'Mapping Read Quality',
    ylabel = 'Mapping Quality (max 60)',
    figsize=(20, 10)
)
QualityFigure = QualityGraph.get_figure()
QualityFigure.savefig(ProjectName + '_QualityGraph.png')
print("QualityFigure Done")

######################################################################### FINAL PRINTOUT #########################################################################

PredictedReference = PredictedReference.replace(".fa","")
ReferenceAccession = PredictedReference.split("_")[0]
ReferenceIsolate = PredictedReference.split("_")[1:]
ReferenceIsolate = " ".join(ReferenceIsolate)

with open(ProjectName + "_SummaryReadout.txt", "w") as f:
    with open("Graphic.txt") as graphicfile:
        f.write(graphicfile.read())
    print("We appreciate your decision to utilize the ASFV assembly genome pipeline, which can be accessed at asfvgenomics.com."
        + "\n" + "Should you have any inquiries concerning this assembly, please feel free to contact us via email at douglas.gladue@usda.gov or genomicsandvaccines@gmail.com." 
        + "\n" 
        + "\n" + "If this service was of use to you please cite: "
        + "\n" + "####TBD#####"
        + "\n" 
        + "\n" + "############################################################################################################"
        + "\n" + WEBSITE 
        + "\n"
        + "\n" 
        + "\n" + "Isolate Name: " + MetaData['isolate']
        + "\n" + "File ID: " + ProjectName
        + "\n" + "Date of Assembly: " + str(datetime.date.today())
        + "\n" + "Predicted Genotype: " + str(Genotype)
        + "\n" + "Predicted Biotype: " + str(Biotype)
        + "\n" + "Closest Match: " + ReferenceIsolate + " (" + str(ReferenceAccession) + ")"
        + "\n" 
        + "\n" + "PROVIDED METADATA"
        + "\n" + "Organism: " + ORGANISM, file = f)
    if MetaData['collection_date'] != None:
        print("Collection Date: " + str(MetaData['collection_date']), file=f)
    if MetaData['country'] != None:
        print("Country: " + str(MetaData['country']), file=f)
    if MetaData['location'] != None:
        print("Location: " + str(MetaData['location']), file=f)
    if MetaData['host'] != None:
        print("Host: " + str(MetaData['host']), file=f)
    if MetaData['tissue_type'] != None:
        print("Tissue Type: " + str(MetaData['tissue_type']), file=f)
    if MetaData['collected_by'] != None:
        print("Collected By: " + str(MetaData['collected_by']), file=f)
    print("\n"
        + "\n" + "############################################################################################################"
        + "\n" + "Simplified pipeline:" 
        + "\n" + "Illumina reads were trimmed and mapped to their closest reference. "
        + "\n" + "A consensus genome was extracted from this mapping and annotated using a currated ASFV CDS database."
        + "\n"
        + "\n" + "Your genome is available in both fasta " + ConsensusGenomeMetaData + " and Genbank " + AnnotatedGenbank + ".gb" + " file formats. "
        + "\n"
        + "\n" + "Illumina reads were mapped back to the consensus to produce the coverage map " + ProjectName + '_CoverageGraph.png'
        + "\n" + "The quality of the mapping can be visualized here: " + ProjectName + '_QualityGraph.png'
        + "\n" + "The raw data that was used to produce these graphs can be found here: " + ProjectName + "_statistics.txt"
        + "\n"
        + "\n" + "Illumina reads were mapped to the reference sequence. Compared to the reference, SNP that were present in over 50% of the reads can be observed here: " + ProjectName + "_SNP.csv"
        + "\n"
        + "\n" + "If you chose to upload this genome, annotations in the required NCBI format can be found here: " + ProjectName + "_FeatureTable.txt"
        + "\n"
        + "\n" + "############################################################################################################"
        + "\n" + "Additional References:"
        + "\n" + "####TBD#####", file = f)


FinalFolder = Make_Directory(ProjectName = ProjectName, Suffix="_TempFiles")
Temp_Folder = Make_Directory(ProjectName = ProjectName, Suffix="_OutPut")

"""Move_File(ProjectName + "_consensus_metadata.fa")
Move_File(ProjectName + "_FeatureTable.txt")
Move_File(ProjectName + "_CoverageGraph.png")
Move_File(ProjectName + "_SummaryReadout.txt")
Move_File(ProjectName + "_statistics.txt")
Move_File(ProjectName + "_SNP.csv")
Move_File(ProjectName + "_QualityGraph.png")
Move_File(ProjectName + "_AnnotationTable.csv")
Move_File(ProjectName + "_Annotated.gb")
Move_File(ProjectName + "_Lane1.html")
Move_File(ProjectName + "_merge_map_consensus_sort.bam")"""
