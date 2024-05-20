import glob
import os
import fnmatch
import pandas as pd
import matplotlib.pyplot as plt
import random
import re
import shutil
import string
import sys
import numpy as np
#
from Bio import SeqIO
from Bio import AlignIO
from Bio import Seq #New Import, worked in IlluminaOnly Pipeline
from Bio.Align import AlignInfo
from Bio.Blast.Applications import NcbiblastxCommandline
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.SeqFeature import SeqFeature, SimpleLocation #New Import, worked in IlluminaOnly Pipeline
import unicodedata #New Import
from statistics import mean

## Combined 3 fasta sequences
def Add_Flanking_Sequence(ProjectName, MiddleFasta_File, x5_PrimeSeqence_File, x3_PrimeSequence_File):
    # Create Final Fasta
    OutFile = ProjectName + "_Consensus_with_flanking_preMap.fa"
    middle_sequence = list(SeqIO.parse(open(MiddleFasta_File),'fasta'))[0].seq
    x5_sequence = list(SeqIO.parse(open(x5_PrimeSeqence_File),'fasta'))[0].seq
    x3_sequence = list(SeqIO.parse(open(x3_PrimeSequence_File),'fasta'))[0].seq
    new_sequences = x5_sequence + middle_sequence + x3_sequence
    with open(OutFile, "w") as f:
        f.write(">" + ProjectName + "\n" + str(new_sequences) + "\n")
    f.close()
    print("***Contig Extender 6/6_Complete_***")
    return OutFile

## Annotate Fasta, create genbank file and NCBI annotation file
def Annotation_Function(ProjectName, Fasta_to_Annotate, Genbank_to_scrape):
    AnnotatedGenbank = ProjectName + "_Annotated"
    Project_GenbankFile = AnnotatedGenbank + ".gb"
    os.system("TheTransporter -t" + Fasta_to_Annotate + " -a " + Genbank_to_scrape + " -o " + AnnotatedGenbank)
    print("***_Annotations Complete_***")
    return Project_GenbankFile

## Convert Annotation to Dataframe - Not Used
def AnnotationToDF(GBKFile, ProjectName = None, Remove_Putative = True):
    GBKFile = SeqIO.read(open(GBKFile, "r"), "genbank")
    if ProjectName == None:
        dfFeatures = pd.DataFrame(columns=['Gene','Sequence', 'Start', 'End', 'Strand'])
    else:
        dfFeatures = pd.DataFrame(columns=['Isolate', 'Gene','Sequence', 'Start', 'End', 'Strand'])
    #
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
        
        # Construct the DataFrame
        if ProjectName == None:
            dfFeatures.loc[len(dfFeatures.index)] = [FeatureName, FeatureSequence, FeatureStart, FeatureEnd, FeatureStrand]
        else:
            dfFeatures.loc[len(dfFeatures.index)] = [ProjectName, FeatureName, FeatureSequence, FeatureStart, FeatureEnd, FeatureStrand]
    # remove any putative genes from dataframe
    if Remove_Putative == True:
        dfFeatures = dfFeatures[dfFeatures['Gene'] != 'putative']
    return dfFeatures

def GenbankDuplicateGeneCheck(Genbank, ProjectName = None):
    DF = AnnotationToDF(Genbank, ProjectName)
    DF = DF[DF.duplicated(['Gene'], keep=False)].sort_values('Gene')
    ImportantDF = DF[~DF['Gene'].str.contains('|'.join(['hypothetical_','MGF_', "ACD_"]))]
    if len(ImportantDF) > 0:
        DuplicateGeneWarning = "Duplicate Genes Found in Genbank Annotation."
    elif len(DF) > 0:
        DuplicateGeneWarning = "Duplicate Genes Found in Genbank Annotation. Only MGFs, Hypotheticals, and ACD genes."
    else:
        DuplicateGeneWarning = None
        DF = None
    return DF, DuplicateGeneWarning

## Apliyes Quality SNP and unmapped regions to referene sequence
def Apply_Variants_to_Consensus(Genome_Fasta, VCF_Index, BED_File, ProjectName, Suffix):
    Consensus_Sequence = ProjectName + "_" + Suffix + ".fa"
    os.system("cat " + Genome_Fasta + " | bcftools consensus " + VCF_Index + " -m " + BED_File + " > " + Consensus_Sequence)
    print("***_Apply Variants to Fasta Complete_***")
    return Consensus_Sequence

## BioType Prediction
def BiotypingWithBlast(queryfile, blastoutputbank, FastaBank, BioTypeReferenceCSV, SaveTempFiles = False):
    GenomeCsv = pd.read_csv(BioTypeReferenceCSV, index_col='Accession')

    Type = True
    # Empty Folder
    my_files_to_Delete = glob.glob(blastoutputbank + '*')
    for file in my_files_to_Delete:
        os.remove(file)
    # Blast Genome Against Fasta Files
    my_Fasta_files = glob.glob(FastaBank + '*.fasta')
    for file in my_Fasta_files:
        filename = file.replace(FastaBank, '').replace('.fasta', '')
        outputlocation = blastoutputbank + filename + '.txt'
        if Type == True:
            blastcmdline = NcbiblastxCommandline(cmd='blastx', query = queryfile, subject= file, max_hsps = 1, out = outputlocation, outfmt = "6 qseqid sseqid evalue pident bitscore length qlen slen qstart qend sstart send qseq sseq", evalue=0.001)
            os.system(str(blastcmdline))

    my_BLAST_output_files = glob.glob(blastoutputbank + '*')
    BM_Pident = pd.DataFrame()
    BM_Ratio = pd.DataFrame()
    BM_PidentRatio = pd.DataFrame()
    for file in my_BLAST_output_files:
        filename = file.replace(blastoutputbank, '').replace('.txt', '')
        prediction = pd.read_csv(file, sep = '\t', names=('qseqid', 'sseqid', 'evalue', 'pident', 'bitscore', 'length', 'qlen', 'slen', 'qstart', 'qend', 'sstart', 'send', 'qseq', 'sseq'))
        for name, pred in prediction.iterrows():
            BM_Pident.loc[pred['sseqid'],[filename]] = pred['pident']/100
            if pred['length']/pred['slen'] > 1:
                Ratio = 1
            else:
                Ratio = pred['length']/pred['slen']
            BM_Ratio.loc[pred['sseqid'],[filename]] = Ratio
            BM_PidentRatio.loc[pred['sseqid'],[filename]] = ((pred['pident']/100)**2)*Ratio

    if SaveTempFiles == True:
        BM_Pident.to_csv("TempPident.csv")
        BM_Ratio.to_csv("TempRatio.csv")
        BM_PidentRatio.to_csv("TempPidentRatio.csv")

    BM_Pident_Mean = BM_Pident.mean(axis=1, numeric_only = True)
    BM_Ratio_Mean = BM_Ratio.mean(axis=1, numeric_only = True)
    BM_PidentRatio_Mean = BM_PidentRatio.mean(axis=1, numeric_only = True)
    BM_Count = BM_PidentRatio.count(axis=1, numeric_only= True)

    Score = (BM_PidentRatio_Mean*(BM_Count**0.5)).sort_values(ascending=False)

    if SaveTempFiles == True:
        Score.to_csv("TempScore.csv")

    B646MatchPercentage = BM_Pident['B646L'][GenomeCsv.loc[Score.index[0]].name]*100
    if B646MatchPercentage >= 95.0:
        if B646MatchPercentage < 98.0:
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
            Statement = " Based on weighted percent identity, the genome exhibited a " + str(BM_PidentRatio_Mean[Score.index[0]]*100 ) + " similarity to its nearest match within the specified Biotype. This percent identity is below our recommended cutoff of 97.5%. Nevertheless, this outcome does not imply that the genome signifies a novel Biotype, however this genome should be analyzed using the complete algorithm described in 'Reclassification of ASFV into 7 Biotypes Using Unsupervised Machine Learning' (Dinhobl et al 2023). Please contact administration for assistance."
        else:
            Statement = ""
        return str(GenomeCsv.loc[Score.index[0]]['Biotype']), Statement
    else:
        Statement = "Warning! B646 match percentage is " + str(B646MatchPercentage) + " which is very low. Ensure proper sequencing and cleaning before submitting."
        BiotypeReturn = ""
        return BiotypeReturn, Statement 

## Combine All Minion Reads
def CombineMinionReads(Directory, ProjectName):
    # Count number of Minion Reads in Directory
    CountofFastqGZ = len(fnmatch.filter(os.listdir(Directory), '*.fastq.gz'))
    # Combine Reads if there are more than one
    if CountofFastqGZ > 1:
        CombinedMinionReadFile = Directory + ProjectName +"_CombinedMinion.fastq.gz"
        os.system("cat " + Directory + "*.fastq.gz > " + CombinedMinionReadFile) 
        print("***_Minion Files Combined Complete_***")
    else:
        CombinedMinionReadFile = Directory + str(fnmatch.filter(os.listdir(Directory), '*.fastq.gz')[0])
        print("***_1 Minon Read Detected, Minion Files Combined Complete_***")
        print(CombinedMinionReadFile)
    return CombinedMinionReadFile, CountofFastqGZ

## Find SNPs, unmapped regions, and extract consensus
def Consensus_Extracter(Genome_Fasta, BAM_File, Max_Depth_Coverage, Suffix, ProjectName):
    #Find all SNPs
    #
    RAW_File = ProjectName + "_" + Suffix + "_calls.vcf"
    os.system("bcftools mpileup -L " + Max_Depth_Coverage + " -Ou -f " + Genome_Fasta + " " + BAM_File + " -d " + Max_Depth_Coverage + " | bcftools call -mv -Ov -o "  + RAW_File)
    print("***_Consensus Extractor 1/3 Complete_***")
    ###############################################################################
    #
    # the first 27 rows can be skipped because they are not formatted properly
    Variant_Table = pd.read_csv(RAW_File, sep="\t", skiprows = 27)
    RAW_Table = pd.read_csv(RAW_File, sep="\t", skiprows = 27)
    
    # Stop Process if there are no Variants
    if len(Variant_Table.index) == 0:
        VCF_index = ""
        CSV_Output = ""
        print("no variants found")
        return VCF_index, CSV_Output
    
    # set file names
    CSV_Output = ProjectName + "_" + Suffix + "_SNP.csv"
    VCF_Output = ProjectName + "_" + Suffix + "_corrected_calls.vcf"
    VCF_index = VCF_Output + ".gz"
    
    # the INFO column has an incredible amount of data that is seperated by ';'
    Variant_Table['INFO'] = Variant_Table['INFO'].str.split(';').fillna(Variant_Table['INFO'])

    # explode the data into multiple rows but keep the original index # on each created row
    Variant_Table=Variant_Table.explode('INFO',ignore_index=False)

    # of the data in info, we only care about 'DP' (total coverage) and DP4 (Ref Read F, Ref Read R, Variant Read F, Variant Read R)
    searchfor = ['DP']
    Variant_Table = Variant_Table[Variant_Table.INFO.str.contains('|'.join(searchfor))]

    # Create another dataframe where DP and DP4 is seperated using '=', name the columns ID and Value
    SplitColumn = Variant_Table['INFO'].str.split(pat = '=', expand = True)
    SplitColumn.columns = ['ID', 'Value']

    # Create DP and DP4 columns, since indeces were preserved, the values for DP and DP4 at the same position will be preserved
    SplitColumn = SplitColumn.pivot(columns = 'ID', values = 'Value')

    # Restore string to integer
    CoverageColumn = SplitColumn['DP'].astype(int)

    # since DP4 is composed of (Ref Read F, Ref Read R, Variant Read F, Variant Read R), we only care about Variant Read F & R
    DP4Column = SplitColumn['DP4'].str.split(',', expand = True).astype(int)

    # The sum of this value gives us the number of reads that contain the variant
    VariantSum = DP4Column[2] + DP4Column[3]

    # The frequency of reads that compare the variants as compared to the total coverage at that position
    Frequency = round(VariantSum / CoverageColumn, 3)
    
    # Clean Up Table
    DropColumns = ['ID', 'FILTER', 'INFO', 'FORMAT', BAM_File]
    Variant_Table = Variant_Table[~Variant_Table.index.duplicated(keep='first')].drop(DropColumns, axis = 1)
    Variant_Table.rename(columns={
        "#CHROM": "Reference_Strain",
        "POS": "Position",
        "REF": "Reference",
        "ALT": "Sequenced_Sample"
        }
    )
    # Create columns with calculated values
    Variant_Table['Count'] = VariantSum
    Variant_Table['Coverage'] = CoverageColumn
    Variant_Table['Frequency'] = Frequency
    
    # Drop Variants less than 50%
    Variant_Table = Variant_Table[Variant_Table.Frequency >= 0.50]
    
    # Drop indices in orginal file that are not in the final SNP Table
    RAW_Table = RAW_Table[RAW_Table.index.isin(Variant_Table.index)]

    # Save Variant File in an User Friendly Format
    Variant_Table.to_csv(CSV_Output, sep=',', index=False)
    
    # Creating A Corrected Variant File for PipeLine
    with open(RAW_File) as myfile:
        first_27_lines = myfile.readlines()[0:27]
    myfile.close()
    with open(VCF_Output,'a') as Corrected_File:
        for line in first_27_lines:
            Corrected_File.write(line)
    Corrected_File.close()
    RAW_Table.to_csv(VCF_Output, sep = '\t', mode='a', index = False, header=True)
    #
    print("***_Consensus Extractor 2/3 Complete_***")
    ##############################################################################
    os.system("bcftools view "+ VCF_Output + " -Oz --write-index -o " + VCF_index)
    print("***_Consensus Extractor 3/3 Complete_***")
    return VCF_index, CSV_Output

## Extract Flanking Minion Reads, Align, Create Consensus, 
def Contig_Extender_for_Polish(ProjectName, MinLengthThreshold, MaxLengthThreshold, MaxSubSample, Threshold, Require_Multiple, MinionMappedFirst10, MinionMappedLast10, First_Output_Fasta, Last_Output_Fasta, First_Output_ALN, Last_Output_ALN):
    # FileOutputs
    x5_Prime_NonPolish_FlankingFile = ProjectName + "_5_prime_non_polish.fa"
    x3_Prime_NonPolish_FlankingFile = ProjectName + "_3_prime_non_polish.fa"
    # Import SAM files as Pandas, first 11 colulmns are mandatory, drop unessary columns
    colnames = ['Qname','Flag','Rname','Pos','MapQ','Cigar','RNext','PNExt','Tlen','Seq','Qual']
    droppedColumns = ['Rname','MapQ','RNext','PNExt','Tlen','Qual']
    First_RAW_Table = pd.read_csv(MinionMappedFirst10, sep="\t", quoting=3, usecols=range(11), header = None, index_col=False)
    Last_RAW_Table = pd.read_csv(MinionMappedLast10, sep="\t", quoting=3, usecols=range(11), header = None, index_col=False)
    First_RAW_Table.columns = colnames
    Last_RAW_Table.columns = colnames
    #
    MinReadNumberForAlignment = 12
    # DropColumns that are not needed
    First_Working_Table = First_RAW_Table.drop(columns=droppedColumns)
    Last_Working_Table = Last_RAW_Table.drop(columns=droppedColumns)
    #
    #########################################################################################################
    #                               Explanation of Regex                                                    #
    # Filter by Regex:                                                                                      #
    #       Regrex for 5' Region of Genome = '^[0-9]{3,4}S'                                                 #
    #              "^" = StartofString "[0-9]" = AnyNumber "{3,4}" Appearing 3 to 4 times "S"               #
    #       Regrex for 3' Region of Genome = '[0-9]{3,4}S$'                                                 #
    #              "[0-9]" = AnyNumber "{3,4}" Appearing 3 to 4 times "S" $= EndofString                    #
    #                                                                                                       #
    #########################################################################################################
    #
    First_Working_Table = First_Working_Table[First_Working_Table.Cigar.str.contains('^[0-9]{3,4}S')].reindex()
    Last_Working_Table = Last_Working_Table[Last_Working_Table.Cigar.str.contains('[0-9]{3,4}S$')].reindex()
    #
    # Determine Flanking Nucleotides by Spliting @ 'S'
    SplitSampleName_First = First_Working_Table.Cigar.str.split(pat = 'S', expand = True)[0]
    #
    ########      Overly complicated method for 3' end because [-1] would not work      #################  
    # write the string mirrored                                                                         #
    Last_Working_Table['Cigar'] = Last_Working_Table['Cigar'].apply(lambda x: x[::-1])                  #
    # replace all letters besides S with #                                                              #
    Last_Working_Table['Cigar'] = Last_Working_Table.Cigar.str.replace('[A-RT-Z]', '#', regex = True)   #
    # Split @ First #                                                                                   #
    Last_Working_Table['Cigar'] = Last_Working_Table.Cigar.str.split(pat = '#', expand = True)[0]       #
    # Correct Mirror                                                                                    #
    Last_Working_Table['Cigar'] = Last_Working_Table['Cigar'].apply(lambda x: x[::-1])                  #
    # Determine Flanking Nucleotides by Spliting @ 'S'                                                  #
    SplitSampleName_Last = Last_Working_Table.Cigar.str.split(pat = 'S', expand = True)[0]              #
    #                                                                                                   #
    #####################################################################################################
    #
    # Create New Column and Convert to Numeric
    First_Working_Table['Flanking'] = SplitSampleName_First
    Last_Working_Table['Flanking'] = SplitSampleName_Last
    First_Working_Table['Flanking'] = pd.to_numeric(First_Working_Table['Flanking'])
    Last_Working_Table['Flanking'] = pd.to_numeric(Last_Working_Table['Flanking'])
    #
    # Calculate the length of Sequence
    First_Working_Table['len'] = First_Working_Table['Seq'].str.len()
    Last_Working_Table['len'] = Last_Working_Table['Seq'].str.len()
    #
    # Extract First N Nucleotides(flanking 5' end) or Last N nucleotides (flanking 3' end)
    First_Working_Table['Seq2'] = First_Working_Table.apply(lambda x: x['Seq'][: x['Flanking'] + 1], 1)
    Last_Working_Table['Seq2'] = Last_Working_Table.apply(lambda x: x['Seq'][x['len'] - x['Flanking'] : x['len'] + 1], 1)
    #
    # Calculate the length of new Sequence
    First_Working_Table['len2'] = First_Working_Table['Seq2'].str.len()
    Last_Working_Table['len2'] = Last_Working_Table['Seq2'].str.len()
    #
    # Filter Reads above certain thereshold
    First_Working_Table = First_Working_Table[First_Working_Table.len2 >= MinLengthThreshold].reindex()
    Last_Working_Table = Last_Working_Table[Last_Working_Table.len2 >= MinLengthThreshold].reindex()
    #
    #
    # Bin and Count the length of the flanking sequences
    bins = range(MinLengthThreshold, MaxLengthThreshold, 500)
    labels = bins[:-1]
    #
    First_Bin_Count = First_Working_Table.groupby(pd.cut(First_Working_Table['len2'], bins=bins, labels=labels)).size().to_frame()
    First_Bin_Count.insert(0, 'LowerBoundary', First_Bin_Count.index)
    First_Bin_Count = First_Bin_Count.rename(columns={0: "Count"})
    # only include bins that contain > MinReadNumberForAlignment
    First_Bin_Count = First_Bin_Count[First_Bin_Count.Count >= MinReadNumberForAlignment].reset_index(drop=True)
    # only perfrom this process if there are enough bins to justify: min length = 1000-1500
    if len(First_Bin_Count) > 1:
        First_Smallest_Size = First_Bin_Count.LowerBoundary[len(First_Bin_Count)- 1]
        First_Largest__Size = First_Smallest_Size + 500
        First_Working_Table = First_Working_Table[First_Working_Table.len2 >= First_Smallest_Size].reindex()
        First_Working_Table = First_Working_Table[First_Working_Table.len2 <= First_Largest__Size].reindex()
        Cancel_5_extend = False
    else:
        Cancel_5_extend = True
        print("5' region will not be extended")
    # Same Process as above but for 3' end
    Last__Bin_Count = Last_Working_Table.groupby(pd.cut(Last_Working_Table['len2'], bins=bins, labels=labels)).size().to_frame()
    Last__Bin_Count.insert(0, 'LowerBoundary', Last__Bin_Count.index)
    Last__Bin_Count = Last__Bin_Count.rename(columns={0: "Count"})
    Last__Bin_Count = Last__Bin_Count[Last__Bin_Count.Count >= MinReadNumberForAlignment].reset_index(drop=True)
    if len(Last__Bin_Count) > 1:
        Last__Smallest_Size = Last__Bin_Count.LowerBoundary[len(Last__Bin_Count)- 1]
        Last__Largest__Size = Last__Smallest_Size + 500
        Last_Working_Table = Last_Working_Table[Last_Working_Table.len2 >= Last__Smallest_Size].reindex()
        Last_Working_Table = Last_Working_Table[Last_Working_Table.len2 <= Last__Largest__Size].reindex()
        Cancel_3_extend = False
    else:
        Cancel_3_extend = True
        print("3' region will not be extended")
    #
    #
    # Extract only a subset of reads to save Computational Time (if needed)
    First_Working_Table_row_Count = First_Working_Table.shape[0]
    Last_Working_Table_row_Count = Last_Working_Table.shape[0]
    if First_Working_Table_row_Count > MaxSubSample:
        First_Working_Table = First_Working_Table.sample(n = MaxSubSample, random_state=1)
    if Last_Working_Table_row_Count > MaxSubSample:
        Last_Working_Table = Last_Working_Table.sample(n = MaxSubSample, random_state=1)
    #
    # Print to CSV (testing purposes only)
    #First_Working_Table.to_csv(First_Output_Fasta[:-2]+"csv", sep=",")
    #Last_Working_Table.to_csv(Last_Output_Fasta[:-2]+"csv", sep=",")
    #
    # Print Flanking Reads to Fasta Files
    if Cancel_5_extend == False:
        First_Name_List = First_Working_Table['Qname'].tolist()
        First_Seq_List = First_Working_Table['Seq2'].tolist()
        with open(First_Output_Fasta, "w") as f:
            for Name, Sequence in zip(First_Name_List, First_Seq_List):
                f.write(">" + Name + "\n" + Sequence + "\n")
        f.close()  
    #
    if Cancel_3_extend == False:
        Last_Name_List = Last_Working_Table['Qname'].tolist()
        Last_Seq_List = Last_Working_Table['Seq2'].tolist()
        with open(Last_Output_Fasta, "w") as f:
            for Name, Sequence in zip(Last_Name_List, Last_Seq_List):
                f.write(">" + Name + "\n" + Sequence + "\n")
        f.close()
    
    #
    print("***_Contig Extender 1/6 Complete_***")
    #
    # Create Alignment with muscle
    if Cancel_5_extend == False:
        os.system("muscle -in " + First_Output_Fasta + " -out " + First_Output_ALN)
    print("***_Contig Extender 2/6_Complete_***")
    if Cancel_3_extend == False:
        os.system("muscle -in " + Last_Output_Fasta + " -out " + Last_Output_ALN)
    print("***_Contig Extender 3/6 Complete_***")
    #
    # Extract Simple Consensus from alignments
    if Cancel_5_extend == False:
        alignment_5prime  = AlignIO.read(open(First_Output_ALN), "fasta")
        consensus_5prime = AlignInfo.SummaryInfo(alignment_5prime).gap_consensus(threshold=Threshold, ambiguous='N', require_multiple=Require_Multiple)
        consensus_5prime = str(consensus_5prime).replace('-','')
    else:
        consensus_5prime = "NNNNN"
    print("***_Contig Extender 4/6 Complete_***")
    #
    if Cancel_3_extend == False:
        alignment_3prime  = AlignIO.read(open(Last_Output_ALN), "fasta")
        consensus_3prime = AlignInfo.SummaryInfo(alignment_3prime).gap_consensus(threshold=Threshold, ambiguous='N', require_multiple=Require_Multiple)
        consensus_3prime = str(consensus_3prime).replace('-','')
    else:
        consensus_3prime = "NNNNN"
    print("***_Contig Extender 5/6 Complete_***")
    #
    # Create Final Fasta
    with open(x5_Prime_NonPolish_FlankingFile, "w") as f:
        f.write(">" + ProjectName + "5_Prime_Flank_nonpolished\n" + str(consensus_5prime) + "\n")
    f.close()
    with open(x3_Prime_NonPolish_FlankingFile, "w") as f:
        f.write(">" + ProjectName + "3_Prime_Flank_nonpolished\n" + str(consensus_3prime) + "\n")
    f.close()
    print("***Contig Extender 6/6_Complete_***")
    return x5_Prime_NonPolish_FlankingFile, x3_Prime_NonPolish_FlankingFile

## Create Depth of Coverage and Mapping Quality Graphs
def Create_Graph(Stats_File, ProjectName):
    plt.close("all")
    StatDataFrame = pd.read_csv(Stats_File, sep='\t', header = None, names = ['Genome', 'Position', 'Coverage', 'Quality'], quoting=3)
    
    # Graph Coverage and Save PNG
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

    # Calculate Average Coverage
    AverageCoverage = round(StatDataFrame.Coverage.mean(), 2)

    # Graph Mapping Quality and Save PNG
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

    print("***_Graphs Complete_***")
    return ProjectName + '_CoverageGraph.png', ProjectName + '_QualityGraph.png', AverageCoverage

## Write Result Summary to Text File
def Create_Results_File(ProjectName, Graphic, Website, isolate, Time, p72_Genotype, p72_Isolate, p72_Accession, p72_PID, p72_Length, p72_Warning, Organism, collection_date, country, location, host, tissue, collected_by, Fasta_Assembly, GenBank_Assembly, Biotype, AverageCoverage, ErrorString, DuplicateGeneWarning):
    with open(ProjectName + "_SummaryReadout.txt", "w") as f:
        """with open(Graphic) as graphicfile:
            f.write(graphicfile.read())
        graphicfile.close()"""
        print("We appreciate your decision to utilize the ASFV assembly genome pipeline. This pipeline was constructed by Edward Spinard, Mark Dinhobl, Douglas Gladue, Jacob Fenster, Nicolos Tesler, Hillary Birtley, Cass Erdelyan, James O'Dwyer, and Anthony Signore."
            + "\n" + "Should you have any inquiries concerning this assembly, please feel free to contact us via email at douglas.gladue@usda.gov." 
            + "\n" 
            + "\n" + "If this service was of use to you please cite: "
            + "\n" + "####TBD#####"
            + "\n" 
            + "\n" + "############################################################################################################"
            + "\n" + Website 
            + "\n"
            + "\n" 
            + "\n" + "Isolate Name: " + isolate
            + "\n" + "File ID: " + ProjectName
            + "\n" + "Date of Assembly: " + Time
            + "\n" + "Average Illumina Read Coverage Across Genome: " + str(AverageCoverage)
            + "\n" + "Predicted Genotype: " + str(p72_Genotype)
            + "\n" + "Closest p72 Match: " + str(p72_Isolate) + " (" + str(p72_Accession) + ")"
            + "\n" + "The B646L(P72) encoded by your genome is " + str(p72_PID) + "%" + " identical to " + str(p72_Length) + ' amino acids encoded by ' + str(p72_Isolate)
            + "\n" + p72_Warning
            + "\n" 
            + "\n" + "Predicted Biotype:" + Biotype
            + "\n" + "PROVIDED METADATA"
            + "\n" + "Organism: " + Organism, file = f)
        if collection_date != None:
            print("Collection Date: " + str(collection_date), file=f)
        if country != None:
            print("Country: " + str(country), file=f)
        if location != None:
            print("Location: " + str(location), file=f)
        if host != None:
            print("Host: " + str(host), file=f)
        if tissue != None:
            print("Tissue Type: " + str(tissue), file=f)
        if collected_by != None:
            print("Collected By: " + str(collected_by), file=f)
        if ErrorString != None:
            print("\n" + "############################################################################################################" + "\n" + "Non-Fatal Errors Detected", file=f)
        if DuplicateGeneWarning != None:
            print(str(DuplicateGeneWarning) + " These Genes can be found at " + ProjectName + "_DuplicateGenes.csv", file=f)
        print("\n"
            + "\n" + "############################################################################################################"
            + "\n" + "Simplified pipeline explanation (Illumina Only):" 
            + "\n" + "Illumina reads were trimmed and mapped to their closest reference."
            + "\n" + "A consensus genome was extracted from this mapping and annotated using a currated ASFV CDS database."
            + "\n" + "Please look over the annotations as gene fusions/splits may not be properly annotated."
            + "\n" + "Without the incorporation of long reads, the detection of large genomic insertions or deletions can be difficult."
            + "\n" + "It is crucial to be mindful that the absence of long reads may impede the accurate identification of these genetic changes."
            + "\n" + "We suggest you visually examine the following read mapping files to confirm any deletions and/or SNPs:"
            + "\n" + ProjectName + "_DenovoScaffold_Mapped_to_" + ProjectName + "_MergedMap_(ReferenceName).sam"
            + "\n" + ProjectName + "_MergedMap_(ReferenceName).bam"
            + "\n"
            + "\n" + "############################################################################################################"
            + "\n" + "Simplified pipeline explanation (Nanopore and Illumina):" 
            + "\n" + "Illumina reads were trimmed and mapped to the ASFV Georgia-2007 (FR682468). Mapped Reads were collected. Unmapped reads were collected and mapped to the genome of Sus scrofa (Sus_scrofa.Sscrofa11.1.dna.toplevel)."
            + "\n" + "Unmapped Illumina reads that remained in pairs were collected and used along with the reads that mapped to ASFV Georgia- 2007 and Nanopore reads to assemble a final contig that was annotated using a currated ASFV CDS database."
            + "\n" + "Please look over the annotations as gene fusions/splits may not be properly annotated."
            + "\n"
            + "\n" + "############################################################################################################"
            + "\n" + "Your genome is available in both fasta " + Fasta_Assembly + " and Genbank " + GenBank_Assembly + " file formats. "
            + "\n"
            + "\n" + "Illumina reads were mapped back to the consensus to produce the coverage map " + ProjectName + '_CoverageGraph.png'
            + "\n" + "The quality of the mapping can be visualized here: " + ProjectName + '_QualityGraph.png'
            + "\n" + "The raw data that was used to produce these graphs can be found here: " + ProjectName + "_statistics.txt"
            + "\n" + "Illumina read mappings can be observed here:" + ProjectName + "_MergedMap_ConsensusGenome.bam"
            + "\n" + "Nanopore read mappings can be observed here:" + ProjectName + "___Minion_Sort.sam"
            + "\n"
            + "\n" + "If you chose to upload this genome, annotations in the required NCBI format can be found here: " + ProjectName + "_FeatureTable.txt"
            + "\n"
            + "\n" + "############################################################################################################"
            + "\n" + "Additional References:"
            + "\n" + "Bankevich, A. et al. SPAdes: a new genome assembly algorithm and its applications to single-cell sequencing. J Comput Biol. 2012 May;19(5):455-77. doi: 10.1089/cmb.2012.0021. Epub 2012 Apr 16. PMID: 22506599; PMCID: PMC3342519."
            + "\n" + "Camacho, C. et al 2009. BLAST+: architecture and applications. BMC Bioinformatics, 10, 421."
            + "\n" + "Cock, P.J. et al. 2009. Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics, 25(11), pp.1422–1423."
            + "\n" + "Dinhobl M, Spinard E, Tesler N, Birtley H, Signore A, Ambagala A, Masembe C, Borca MV, Gladue DP. Reclassification of ASFV into 7 Biotypes Using Unsupervised Machine Learning. Viruses. 2023 Dec 30;16(1):67. doi: 10.3390/v16010067. PMID: 38257767; PMCID: PMC10819123."
            + "\n" + "Dinhobl M, Spinard E, Birtley H, Tesler N, Borca MV, Gladue DP. African swine fever virus P72 genotyping tool. Microbiol Resour Announc. 2024 Feb 15;13(2):e0089123. doi: 10.1128/mra.00891-23. Epub 2024 Jan 8. PMID: 38189309; PMCID: PMC10868209."
            + "\n" + "Edgar RC. MUSCLE: multiple sequence alignment with high accuracy and high throughput. Nucleic Acids Res. 2004 Mar 19;32(5):1792-7. doi: 10.1093/nar/gkh340. PMID: 15034147; PMCID: PMC390337."
            + "\n" + "Harris, C.R. et al. Array programming with NumPy. Nature 585, 357–362 (2020). DOI: 10.1038/s41586-020-2649-2."
            + "\n" + "Hunter, J. D. Matplotlib: A 2D Graphics Environment, in Computing in Science & Engineering, vol. 9, no. 3, pp. 90-95, May-June 2007, doi: 10.1109/MCSE.2007.55."
            + "\n" + "Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics, 34:3094-3100. doi:10.1093/bioinformatics/bty191"
            + "\n" + "Li, H. (2021). New strategies to improve minimap2 alignment accuracy. Bioinformatics, 37:4572-4574. doi:10.1093/bioinformatics/btab705"
            + "\n" + "McKinney, W. et al. 2010. Data structures for statistical computing in python. In Proceedings of the 9th Python in Science Conference. pp. 51–56."
            + "\n" + "Danecek, P. et al. Twelve years of SAMtools and BCFtools, GigaScience, Volume 10, Issue 2, February 2021, giab008, https://doi.org/10.1093/gigascience/giab008"
            + "\n" + "Quinlan AR and Hall IM, 2010. BEDTools: a flexible suite of utilities for comparing genomic features. Bioinformatics. 26, 6, pp. 841–842."
            + "\n" + "Chen, S. 2023. Ultrafast one-pass FASTQ data preprocessing, quality control, and deduplication using fastp. iMeta 2: e107. https://doi.org/10.1002/imt2.107"
            + "\n" + "Chen, S. et al. ; fastp: an ultra-fast all-in-one FASTQ preprocessor, Bioinformatics, Volume 34, Issue 17, 1 September 2018, Pages i884–i890, https://doi.org/10.1093/bioinformatics/bty560"
            + "\n" + "Spinard E, Dinhobl M, Tesler N, Birtley H, Signore AV, Ambagala A, Masembe C, Borca MV, Gladue DP. A Re-Evaluation of African Swine Fever Genotypes Based on p72 Sequences Reveals the Existence of Only Six Distinct p72 Groups. Viruses. 2023 Nov 11;15(11):2246. doi: 10.3390/v15112246. PMID: 38005923; PMCID: PMC10675559."
            + "\n" + "Vasimuddin, Md. et al. Efficient Architecture-Aware Acceleration of BWA-MEM for Multicore Systems. IEEE Parallel and Distributed Processing Symposium (IPDPS), 2019. 10.1109/IPDPS.2019.00041", file = f)
    f.close()
    return ProjectName + "_SummaryReadout.txt"

## Delete Files
def Delete_File(DeleteFile):
    if DeleteFile != None:
        if os.path.exists(DeleteFile):
            os.remove(DeleteFile)
        else:
            print(str(DeleteFile + "does not exist"))
    else:
        print("empty variable") 

## Delete Folders
def Delete_Folder(DeleteFolder):
    if DeleteFolder != None:
        try:
            shutil.rmtree(DeleteFolder)
        except OSError as e:
            print("Error: %s - %s." % (e.filename, e.strerror))
    else:
        print("empty variable")

## End Script
def End_Script(Total_Contig_Length):
    if Total_Contig_Length < 150000:
        return True
    else:
        return False

## Extract from BAM File
def ExtractFromMerged(Input, Reference, ProjectName):
    
    # Explanation:
    #   -b  =   Output in the BAM format. 
    #   -f  =   Only output alignments with all bits set in INT present in the FLAG field
    #   -F  =   Do not output alignments with any bits set in INT present in the FLAG field
    #   0x2 = 	PROPER_PAIR	each segment properly aligned according to the aligner
    os.system("samtools view -b -f 0x2 " + Input + " | samtools sort -n -o " + ProjectName + "_mappedReads_" + Reference + "_sort.bam")
    print("***_Reads Extracted Complete 1/2_***")

    #extract as fastq-paired Reads
    ExtractedRead01 = ProjectName + "_" + Reference + "_TrimmedExtracted1.fastq.gz"
    ExtractedRead02 = ProjectName + "_" + Reference + "_TrimmedExtracted2.fastq.gz"
    
    os.system("samtools fastq -@ 8 " + ProjectName + "_mappedReads_" + Reference + "_sort.bam -1 " + ExtractedRead01 + " -2 " + ExtractedRead02 + " -0 /dev/null -s /dev/null -n")
    
    print("***_Reads Extracted Complete 2/2_***")
    return ExtractedRead01, ExtractedRead02

## Extract Reads that Did NOT Map
def ExtractFromMerged_UnMappedReads(Input, Reference, ProjectName):
    
    # Explanation:
    #   -b  =   Output in the BAM format. 
    #   -f  =   Only output alignments with all bits set in INT present in the FLAG field
    #   -F  =   Do not output alignments with any bits set in INT present in the FLAG field
    #   0x2 = 	PROPER_PAIR	each segment properly aligned according to the aligner
    #   0x4 =	UNMAP, segment unmapped
    os.system("samtools view -b -f 0x4 " + Input + " | samtools sort -n -o " + ProjectName + "_mappedReads_" + Reference + "_sort.bam")
    print("***_Reads Extracted Complete 1/2_***")

    #extract as fastq-paired Reads
    ExtractedRead01 = ProjectName + "_" + Reference + "_TrimmedExtracted1.fastq.gz"
    ExtractedRead02 = ProjectName + "_" + Reference + "_TrimmedExtracted2.fastq.gz"
    
    os.system("samtools fastq -@ 8 " + ProjectName + "_mappedReads_" + Reference + "_sort.bam -1 " + ExtractedRead01 + " -2 " + ExtractedRead02 + " -0 /dev/null -s /dev/null -n")
    
    print("***_Reads Extracted Complete 2/2_***")
    return ExtractedRead01, ExtractedRead02

## Extract Largest Contig
def ExtractLargestContig(Assembled_Scaffold_File, ProjectName):
    AssembledContigs = []
    GenomeFileTemp = ProjectName + "_temp_assembly.fa"

    for seq_record in SeqIO.parse(Assembled_Scaffold_File, "fasta"):
        AssembledContigs.append([seq_record.id, str(seq_record.seq), len(seq_record)])
    
    # Sort Scaffolds by length
    AssembledContigs.sort(key=lambda x: x[2])

    # Extract Longest ID and Sequence
    LongestFasterHeader = AssembledContigs[-1][0]
    LongestSequence= AssembledContigs[-1][1]
    del AssembledContigs

    # Variable to be used later
    Total_Contig_Length = int(LongestFasterHeader.split("_")[3])
    

    # Write Data to file
    with open(GenomeFileTemp, "w") as f:
        f.write(">" + ProjectName + "\n" + str(LongestSequence) + "\n")
    f.close()

    print("***_Largest Contig Extracted Complete_***")
    return GenomeFileTemp, Total_Contig_Length

## Process BAM file, Create Stats file
def Generate_Mapping_Stats(GenomeFile, ProjectName):
    SortBamOutput = ProjectName + "_merge_map_consensus_sort.bam"
    ConsensusGenome_Stats_File = ProjectName + "_statistics.txt"
    os.system("samtools view -b -f 0x2 " + GenomeFile + " | samtools sort -o " + SortBamOutput) #-n flag will cause mpileup to fail
    print("***_Mapping Stats 1/2 Complete_***")
    os.system("samtools mpileup " + ProjectName + "_merge_map_consensus_sort.bam" + " -s -a -d 0 | cut -f3,5,6 --complement > " + ConsensusGenome_Stats_File)
    print("***_Mapping Stats 2/2 Complete_***")
    return ConsensusGenome_Stats_File

## Index Genome
def Genome_Index(Genome_To_Index):
    os.system("bwa-mem2 index " + Genome_To_Index)
    print("***_Genome Index_***")

## Map Illumina Reads to Reference
def MapToReference(Trimmed_Read01, Trimmed_Read02, Lane, Threads, ReferenceDirectory, Reference, ProjectName, Stringent = False):
    if Trimmed_Read01 == None or Trimmed_Read02 == None:
        OutputLocation = ""
        return OutputLocation
    
    if Stringent == True:
        K_Parameter = " -k 45 "
    else:
        K_Parameter = " "

    ReferenceLocation = ReferenceDirectory + Reference
    OutputLocation = ProjectName + "_mapped_" + Reference.replace(".fa","") + "_Lane" + str(Lane) + ".bam"

    os.system("bwa-mem2 mem -t " + Threads + K_Parameter + ReferenceLocation + " " + Trimmed_Read01 + " " + Trimmed_Read02 + " | samtools sort -o " + OutputLocation + " " + "-")
    
    print("***_Illumina Mapping Complete_***")
    return OutputLocation

## MiniMap Minion Reads to Largest Contig and Extract Reads (Note: "MiniMapToReference:samtools view" differs from function "ExtractFromMerged: samtools view")
def MiniMapToReference(GenomeFileTemp, CombinedMinionReadFile, ProjectName, x5_Start_NT, x5_End_NT, x3_Start_NT, x3_End_NT):
    MinionMappedTempSAM = ProjectName + "_Minion_Mapped.sam"
    MinionMappedTempBAM = ProjectName + "_Minion_Mapped.bam"
    MinionMappedTempSOR = ProjectName + "___Minion_Sort.bam"
    MinionMappedFirst10 = ProjectName + "_first_10.sam"
    MinionMappedLast10 = ProjectName + "_last_10.sam"
    os.system("minimap2 -ax map-ont " + GenomeFileTemp + " " + CombinedMinionReadFile + " -o " + MinionMappedTempSAM)
    print("***_MiniMap 1/6 Complete_***")
    os.system("samtools view -b -F 4 " + MinionMappedTempSAM + " > " + MinionMappedTempBAM)
    print("***_MiniMap 2/6 Complete_***")
    os.system("samtools sort " + MinionMappedTempBAM + " -o " + MinionMappedTempSOR)
    print("***_MiniMap 3/6 Complete_***")
    os.system("samtools index " + MinionMappedTempSOR)
    print("***_MiniMap 4/6 Complete_***")
    os.system("samtools view " + MinionMappedTempSOR + " " + ProjectName + ":" + x5_Start_NT + "-" + x5_End_NT + " > " + MinionMappedFirst10)
    print("***_MiniMap 5/6 Complete_***")
    os.system("samtools view " + MinionMappedTempSOR + " " + ProjectName + ":" + x3_Start_NT + "-" + x3_End_NT + " > " + MinionMappedLast10)
    print("***_MiniMap 6/6 Complete_***")
    return [MinionMappedTempSAM, MinionMappedTempBAM, MinionMappedTempSOR, MinionMappedFirst10, MinionMappedLast10]

## Combine All BAM Mappings
def MergeMaps(MapL1, MapL2, MapL3, MapL4, ReferenceName, ProjectName):
    if str(MapL2) + str(MapL3) + str(MapL4) == "":
        CombinedMap = MapL1
    else:
        CombinedMap = ProjectName + "_MergedMap_" + ReferenceName + ".bam"
        Command = "samtools merge " + CombinedMap + " " + str(MapL1) + " " +  str(MapL2) + " " + str(MapL3) + " " + str(MapL4)
        os.system(Command.replace("  ", " ").replace("  ", " ").replace("  ", " "))
    print("***_Merged Mappings Complete_***")
    return CombinedMap

## Add Metadata to final Assembly
def MetaDataString(input_fasta_file, isolate, organism, collection_date, country, location, host, tissue, collected_by, Genotype, Biotype, note):
    #
    FastaFile = list(SeqIO.parse(input_fasta_file, format = 'fasta'))
    ConsensusGenomeMetaData = input_fasta_file.replace('.fa', '_metadata.fa')
    ConsensusGenomeMetaData = input_fasta_file.replace('preMeta', '')
    string = str()
    #
    if isolate != None:
        string += isolate + " [Isolate=" + isolate + "]"
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
    if Genotype != None:
        string += " [Genotype=" + str(Genotype) + "]"
    if Biotype != None:
        string += " [Biotype=" + Biotype + "]"
    if note != None:
        string += " [Note=" + note + "]"
    #
    with open(ConsensusGenomeMetaData, "w") as f:
        print(">" + string, file = f)
        print(str(FastaFile[0].seq), file = f)
        f.close()

    print("***_Metadata Header added to Fasta Complete_***")
    return ConsensusGenomeMetaData

#TODO: NEW FUNCTION
def MetaDataRead(NewMetaData):    
    NewMetaData = pd.read_csv(NewMetaData)
    NewMetaData = NewMetaData.where(pd.notnull(NewMetaData), None)
    NewMetaData = NewMetaData.replace(np.nan, None)
    return NewMetaData

## Update MetadataFile
def MetaDataUpdate(MetaDataFile, NewMetaData):

    MetaDataDF = pd.read_csv(MetaDataFile)
    MetaDataDF = MetaDataDF.where(pd.notnull(MetaDataDF), None)

    NewMetaData['Errors'] = "Failure. Uknown: Run did not Complete."
    MetaDataToStore = pd.concat([MetaDataDF,NewMetaData.to_frame().transpose()],axis=0)
    MetaDataToStore.to_csv(MetaDataFile, index=False)

    print("***_Metadata Upated_***")

def ErrorUpdate(MetaDataFile, Error, ProjectName):
    
    MetaDataDF = pd.read_csv(MetaDataFile, index_col="Project_ID (internal)")
    MetaDataDF = MetaDataDF.where(pd.notnull(MetaDataDF), None)
    MetaDataDF.at[ProjectName, "Errors"] = Error
    MetaDataDF.to_csv(MetaDataFile)
    print(Error)

## Create a new directory
def Make_Directory(ProjectName, SecondPart = "_BLASTOUTPUTBANK/"):
    NewDirectory = ProjectName + SecondPart
    if not os.path.exists(NewDirectory):
        os.makedirs(NewDirectory)
        return NewDirectory
    else:
        return NewDirectory

## Move Files
def Move_File(ProjectName, List_of_Files, Title = "Project"):
    NewDirectory = Title + "_" + str(ProjectName) +"/"
    if not os.path.exists(NewDirectory):
        os.makedirs(NewDirectory)
    else:
        pass
    for file in List_of_Files:
        try:
            if file != None:
                if os.path.exists(file):
                    shutil.move(file, NewDirectory + file)
            else:
                print("empty variable")
        except:
            print("Failed to move " + file)

## Predict p72 Genotype, Updated to return query start and query end
def p72finder_new(queryfile, subjectfile, ProjectName):
    #
    GenomeColumn = 'Genome'
    SequenceColumn = 'Sequence'
    IsolateColumn = 'Isolate'
    GenotypeColumn = 'GenotypeNumber'
    outputname_fasta = ProjectName + '_P72.fasta'
    outputlocation = ProjectName + '_RAW_p72_Genotype_blastxout.txt'
    #
    subject_excel = pd.read_csv(subjectfile, usecols=[GenomeColumn,SequenceColumn,IsolateColumn,GenotypeColumn])

    with open(outputname_fasta, "w") as f:
        pass

    for index, row in subject_excel.iterrows():      
        with open(outputname_fasta, "a") as f:
            print(">" + row[GenomeColumn] + "\n" + str(row[SequenceColumn]).replace("[", "").replace("]", "").replace("'", ""), file = f)

    # Identify if Nucleotide or Amino Acid Sequence
    #for record in Bio.SeqIO.parse(queryfile, "fasta"): old line replacedd
    for record in SeqIO.parse(queryfile, "fasta"):
        ok = r'ACTGUWSMKRYBDHVN*'
        Type = all(c in ok for c in str(record.seq))
    
    if Type == True:
        from Bio.Blast.Applications import NcbiblastxCommandline
        blastcmdline = NcbiblastxCommandline(cmd='blastx', query = queryfile, subject= outputname_fasta, max_hsps = 1, max_target_seqs = 43, out = outputlocation, outfmt = "6 qseqid sseqid evalue pident bitscore length qlen slen qstart qend sstart send qseq sseq", evalue=0.001)
        os.system(str(blastcmdline))
    if Type == False:
        from Bio.Blast.Applications import NcbiblastpCommandline
        blastcmdline = NcbiblastpCommandline(cmd='blastp', query = queryfile, subject= outputname_fasta, max_hsps = 1, max_target_seqs = 43, out = outputlocation, outfmt = "6 qseqid sseqid evalue pident bitscore length qlen slen qstart qend sstart send qseq sseq", evalue=0.001)
        os.system(str(blastcmdline))

    # Read Blast Results
    prediction = pd.read_csv(outputlocation, sep = '\t', names=('qseqid', 'sseqid', 'evalue', 'pident', 'bitscore', 'length', 'qlen', 'slen', 'qstart', 'qend', 'sstart', 'send', 'qseq', 'sseq')).drop_duplicates(["qseq"])
    result = prediction.merge(subject_excel, left_on='sseqid', right_on=GenomeColumn)
    
    MaxBitscore = result['bitscore'].max()
    MaxPident = result[result['bitscore'] == MaxBitscore]['pident'].max()

    result = result[(result['bitscore'] == MaxBitscore) & (result['pident'] == MaxPident)]
    Predicted_Genotype = str(result.iloc[0][GenotypeColumn])
    Predicted_Genotype_PID = str(result.iloc[0]['pident'])
    Predicted_Genotype_MatchLength = str(result.iloc[0]['length'])
    Predicted_Genotype_Match_Isolate = str(result.iloc[0][IsolateColumn])
    Predicted_Genotype_Match_Accession = str(result.iloc[0][GenomeColumn])
    Predicted_Genotype_qstart = str(result.iloc[0]['qstart'])
    Predicted_Genotype_qend = str(result.iloc[0]['qend'])
    
    print('Predicted Genotype: ' + Predicted_Genotype + "\n" + "The B646L(P72) encoded by your genome is " + Predicted_Genotype_PID + "%" + " identical to the " + Predicted_Genotype_MatchLength + ' amino acids in ' + Predicted_Genotype_Match_Isolate + ' (' + Predicted_Genotype_Match_Accession + ')')
    
    if result.iloc[0]['length'] <= 600:
        WarningString = "WARNING!!! No full length B646L found in the submitted sequence. A manual double-check of raw blast output is recommended, and genotyping results are highly likely to be illegitimate."
    elif result.iloc[0]['length'] <= 635:
        WarningString = "Warning - No full length B646L found in submitted sequence. A manual double-check of the raw blast output is recommended, there are potentially multiple genotypes overlapping."
    else:
        WarningString = ""
    print("***_Genotyping Complete_***")
    return [Predicted_Genotype, Predicted_Genotype_PID, Predicted_Genotype_MatchLength, Predicted_Genotype_Match_Isolate, Predicted_Genotype_Match_Accession, WarningString, outputlocation, Predicted_Genotype_qstart, Predicted_Genotype_qend]

## Reverse Complements Genome if p72 is in the incorrect orientation
def QStart_Qend_NewFasta(queryfile, qstart, qend):
    if float(qstart) < float(qend):
        with open(queryfile) as fasta_file:
            for seq_record in SeqIO.parse(fasta_file, 'fasta'):
                    desc = seq_record.description
                    seq = seq_record.seq.reverse_complement() 
            fasta_file.close()
                    
        with open(queryfile, "w") as f:
            print(">" + desc + "\n" + seq, file = f)
            f.close()       

def p72finder(queryfile, subjectfile, ProjectName):
    #
    GenomeColumn = 'Genome'
    SequenceColumn = 'Sequence'
    IsolateColumn = 'Isolate'
    GenotypeColumn = 'GenotypeNumber'
    outputname_fasta = ProjectName + '_P72.fasta'
    outputlocation = ProjectName + '_RAW_p72_Genotype_blastxout.txt'
    #
    subject_excel = pd.read_csv(subjectfile, usecols=[GenomeColumn,SequenceColumn,IsolateColumn,GenotypeColumn])

    with open(outputname_fasta, "w") as f:
        pass

    for index, row in subject_excel.iterrows():      
        with open(outputname_fasta, "a") as f:
            print(">" + row[GenomeColumn] + "\n" + str(row[SequenceColumn]).replace("[", "").replace("]", "").replace("'", ""), file = f)

    # Identify if Nucleotide or Amino Acid Sequence
    #for record in Bio.SeqIO.parse(queryfile, "fasta"): old line replacedd
    for record in SeqIO.parse(queryfile, "fasta"):
        ok = r'ACTGUWSMKRYBDHVN*'
        Type = all(c in ok for c in str(record.seq))
    
    if Type == True:
        from Bio.Blast.Applications import NcbiblastxCommandline
        blastcmdline = NcbiblastxCommandline(cmd='blastx', query = queryfile, subject= outputname_fasta, max_hsps = 1, max_target_seqs = 43, out = outputlocation, outfmt = "6 qseqid sseqid evalue pident bitscore length qlen slen qstart qend sstart send qseq sseq", evalue=0.001)
        os.system(str(blastcmdline))
    if Type == False:
        from Bio.Blast.Applications import NcbiblastpCommandline
        blastcmdline = NcbiblastpCommandline(cmd='blastp', query = queryfile, subject= outputname_fasta, max_hsps = 1, max_target_seqs = 43, out = outputlocation, outfmt = "6 qseqid sseqid evalue pident bitscore length qlen slen qstart qend sstart send qseq sseq", evalue=0.001)
        os.system(str(blastcmdline))

    # Read Blast Results
    prediction = pd.read_csv(outputlocation, sep = '\t', names=('qseqid', 'sseqid', 'evalue', 'pident', 'bitscore', 'length', 'qlen', 'slen', 'qstart', 'qend', 'sstart', 'send', 'qseq', 'sseq')).drop_duplicates(["qseq"])
    result = prediction.merge(subject_excel, left_on='sseqid', right_on=GenomeColumn)
    
    MaxBitscore = result['bitscore'].max()
    MaxPident = result[result['bitscore'] == MaxBitscore]['pident'].max()

    result = result[(result['bitscore'] == MaxBitscore) & (result['pident'] == MaxPident)]
    Predicted_Genotype = str(result.iloc[0][GenotypeColumn])
    Predicted_Genotype_PID = str(result.iloc[0]['pident'])
    Predicted_Genotype_MatchLength = str(result.iloc[0]['length'])
    Predicted_Genotype_Match_Isolate = str(result.iloc[0][IsolateColumn])
    Predicted_Genotype_Match_Accession = str(result.iloc[0][GenomeColumn])
    
    print('Predicted Genotype: ' + Predicted_Genotype + "\n" + "The B646L(P72) encoded by your genome is " + Predicted_Genotype_PID + "%" + " identical to the " + Predicted_Genotype_MatchLength + ' amino acids in ' + Predicted_Genotype_Match_Isolate + ' (' + Predicted_Genotype_Match_Accession + ')')
    
    if result.iloc[0]['length'] <= 600:
        WarningString = "WARNING!!! No full length B646L found in the submitted sequence. A manual double-check of raw blast output is recommended, and genotyping results are highly likely to be illegitimate."
    elif result.iloc[0]['length'] <= 635:
        WarningString = "Warning - No full length B646L found in submitted sequence. A manual double-check of the raw blast output is recommended, there are potentially multiple genotypes overlapping."
    else:
        WarningString = ""
    print("***_Genotyping Complete_***")
    return [Predicted_Genotype, Predicted_Genotype_PID, Predicted_Genotype_MatchLength, Predicted_Genotype_Match_Isolate, Predicted_Genotype_Match_Accession, WarningString, outputlocation]

## Generate Random String and Time
def rand_pass(size): 
    generate_pass = ''.join([random.choice(string.ascii_lowercase + string.digits) for n in range(size)]) 
    return generate_pass

## Rename Files
def Rename_File(OldFile, NewFile):
    os.rename(OldFile, NewFile)

## Trim Ambigous Nucleotides at the start and end of a fasta file
def Trim_N_at_Start_End(input_fasta, ProjectName):
    output_fasta = ProjectName + "_preMeta_consensus.fa"

    record = SeqIO.read(input_fasta,'fasta')
    record_id = str(record.id)
    record_sequence = str(record.seq)

    Start = re.search(r'[^N]', record_sequence).start()
    End = re.search(r'[^N]', record_sequence[::-1]).start()

    if End != 0:
        record_sequence = record_sequence[Start:-End]
    else:
        record_sequence = record_sequence[Start:]
        
    with open(output_fasta, "w") as f:
        print(">" + record_id + "\n" + record_sequence, file = f)
    f.close()
    return output_fasta

## Trim Illumina Reads (Do not Trim Minion)
def TrimReads(Threads, ProjectName, Read01, Read02, Lane = 1, Trim_5_Prim = "25", Trim_3_Prim = "5", Ambig = "2", MinLength = "50", QualityPhred = "20", TrimmedRead01Output = None, TrimmedRead02Output = None):
    
    if Read01 == None or Read02 == None:
        Trimmed_Read01 = None
        Trimmed_Read02 = None
        HTML_File = None
        return Trimmed_Read01, Trimmed_Read02, HTML_File
    
    #Setting Output Names
    if TrimmedRead01Output == None:
        Trimmed_Read01 = ProjectName + "_Lane" + str(Lane) + "_Trimmed_R1.fastq.gz"
    else:
        Trimmed_Read01 = TrimmedRead01Output
    if TrimmedRead02Output == None:
        Trimmed_Read02 = ProjectName + "_Lane" + str(Lane) + "_Trimmed_R2.fastq.gz"
    else:
        Trimmed_Read02 = TrimmedRead02Output

    HTML_File = ProjectName + "_Lane" + str(Lane) + ".html"
    os.system ("fastp -V --detect_adapter_for_pe -i " + Read01 + " -I " + Read02 + " -o " + Trimmed_Read01 + " -O " + Trimmed_Read02 + " -f " + Trim_5_Prim + " -F " + Trim_5_Prim + " -t " + Trim_3_Prim + " -T " + Trim_3_Prim + " -n " + Ambig + " -l " + MinLength + " -w " + Threads + " -q " + QualityPhred + " --html " + HTML_File)
    
    print("***_Trim Reads Complete_***")
    return Trimmed_Read01, Trimmed_Read02, HTML_File

## Find UnMapped regions
def UnMappedRegions(BAM_File, Suffix, ProjectName):
    UnMappedRegions_Bed_File = ProjectName + "_unmappedRegions_" + Suffix + ".bed"
    os.system("bedtools genomecov -ibam " + BAM_File + " -bga | grep -w 0$ > " + UnMappedRegions_Bed_File)
    print("***_UnMapped Regions Complete_***")
    return UnMappedRegions_Bed_File

## Find Regions with Poor Quality Mapping

def FindPoorQualityRegions(statistics_file, Quality=59, Streak=100):
        
    """
    Find low quality regions of a _statistics file, such as one produced by samtools mpileup. 

    Parameters
    ---
    statistics_file : The path to a .txt _statistics file, or a pd.dataframe constructed as one. Quality in the df.dataframe can be either in Phred Quality scores or their translated versions.
    
    >>> Example
    "21g5i9hj_2024_01_30_statistics.txt"
    
    >>> Example
        Genome          Position  Coverage  Quality
    0   Example_Name    1          4        0.0
    1   Example_Name    2          4        0.0
    2   Example_Name    3          4        0.0
    
    Quality: An integer or float. Maxiumum Average (by base pair) Quality Desired for Region
    
    Streak: A positive integer. Length of positions with low quality in a row before being declared as a low quality.

    Returns
    ---
    A list of two-length lists detailing the boundaries of poorly sequenced regions.

    Example
    ---
    >>> df = FindPoorQualityRegions("21g5i9hj_2024_01_30_statistics.txt")
    """
    if type(statistics_file) is pd.DataFrame:
        try: 
            statistics_file['Quality'] = statistics_file['Quality'].apply(lambda x: round(mean([ord(i) for i in x]), 1) - 33)
        except:
            try:
                statistics_file['Mapping_Quality'] = statistics_file['Mapping_Quality'].apply(lambda x: round(mean([ord(i) for i in x]), 1) - 33)
            except:
                pass
    else:
        try:
            statistics_file = pd.read_csv(statistics_file, sep='\t', header = None, names = ['Genome', 'Position', 'Coverage', 'Quality'], quoting=3)
            statistics_file['Quality'] = statistics_file['Quality'].apply(lambda x: round(mean([ord(i) for i in x]), 1) - 33)
        except:
            statistics_file = pd.read_csv(statistics_file, sep='\t', header = None, names = ['Genome', 'Position', 'Coverage', 'Mapping_Quality', 'Quality_Cigar'], quoting=3)
            statistics_file['Mapping_Quality'] = statistics_file['Mapping_Quality'].apply(lambda x: round(mean([ord(i) for i in x]), 1) - 33)

    def getStreaks(L,N,minSize):
        breaks = [i for i,(a,b) in enumerate(zip(L,L[1:]),1) if (a<=N)!=(b<=N)]
        return [[s,e-1] for s,e in zip([0]+breaks,breaks+[len(L)])
                            if e-s>=minSize and L[s]<=N]
    
    try:
        return getStreaks(statistics_file['Quality'], Quality, Streak), statistics_file
    except:
        return getStreaks(statistics_file['Mapping_Quality'], Quality, Streak), statistics_file

def FixGenbankMakeFeatureTable(ProjectName, GBKFile, ConcensusFasta, statistics_file = None):

    #Identify Poor Quality Regions
    def FindPoorQualityRegions(statistics_file, Quality=59, Streak=100):
        
        """
        Find low quality regions of a _statistics file, such as one produced by samtools mpileup. 

        Parameters
        ---
        statistics_file : The path to a .txt _statistics file, or a pd.dataframe constructed as one. Quality in the df.dataframe can be either in Phred Quality scores or their translated versions.
        
        >>> Example
        "21g5i9hj_2024_01_30_statistics.txt"
        
        >>> Example
            Genome          Position  Coverage  Quality
        0   Example_Name    1          4        0.0
        1   Example_Name    2          4        0.0
        2   Example_Name    3          4        0.0
        
        Quality: An integer or float. Maxiumum Average (by base pair) Quality Desired for Region
        
        Streak: A positive integer. Length of positions with low quality in a row before being declared as a low quality.

        Returns
        ---
        A list of two-length lists detailing the boundaries of poorly sequenced regions.

        Example
        ---
        >>> df = FindPoorQualityRegions("21g5i9hj_2024_01_30_statistics.txt")
        """
        if type(statistics_file) is pd.DataFrame:
            try: 
                statistics_file['Quality'] = statistics_file['Quality'].apply(lambda x: round(mean([ord(i) for i in x]), 1) - 33)
            except:
                try:
                    statistics_file['Mapping_Quality'] = statistics_file['Mapping_Quality'].apply(lambda x: round(mean([ord(i) for i in x]), 1) - 33)
                except:
                    pass
        else:
            try:
                statistics_file = pd.read_csv(statistics_file, sep='\t', header = None, names = ['Genome', 'Position', 'Coverage', 'Quality'], quoting=3)
                statistics_file['Quality'] = statistics_file['Quality'].apply(lambda x: round(mean([ord(i) for i in x]), 1) - 33)
            except:
                statistics_file = pd.read_csv(statistics_file, sep='\t', header = None, names = ['Genome', 'Position', 'Coverage', 'Mapping_Quality', 'Quality_Cigar'], quoting=3)
                statistics_file['Mapping_Quality'] = statistics_file['Mapping_Quality'].apply(lambda x: round(mean([ord(i) for i in x]), 1) - 33)

        def getStreaks(L,N,minSize):
            breaks = [i for i,(a,b) in enumerate(zip(L,L[1:]),1) if (a<=N)!=(b<=N)]
            return [[s,e-1] for s,e in zip([0]+breaks,breaks+[len(L)])
                                if e-s>=minSize and L[s]<=N]
        
        try:
            return getStreaks(statistics_file['Quality'], Quality, Streak), statistics_file
        except:
            return getStreaks(statistics_file['Mapping_Quality'], Quality, Streak), statistics_file

    if statistics_file != None:
        PoorQualityRegions, StatDF_Check = FindPoorQualityRegions(statistics_file)

    #Step1: Add Concensus Data to Record

    for record in SeqIO.parse(ConcensusFasta, "fasta"):
        record.annotations = { "molecule_type" : "DNA" }
        
    with open(ProjectName + "_FeatureTable.txt", "w") as f:
        f.write(">Feature " + record.id + "\n")
        f.close()

    #Step2: Clean up Data Entries
    GBKFile = SeqIO.read(open(GBKFile, "r"), "genbank")
    
    #Add Repeat Regions to GB
    if statistics_file != None:
        
        RepeatRegionFound_5Prime = False #Setting to False if they are not found
        RepeatRegionFound_3Prime = False #Setting to False if they are not found
        
        def overlap(start1, end1, start2, end2):
            """how much does the range (start1, end1) overlap with (start2, end2)"""
            return max(max((end2-start1), 0) - max((end2-end1), 0) - max((start2-start1), 0), 0)
        
        if len(PoorQualityRegions) >= 1:
            if overlap(0, 100, PoorQualityRegions[0][0], PoorQualityRegions[0][1]) >= 1:
                feature = SeqFeature(SimpleLocation(PoorQualityRegions[0][0], PoorQualityRegions[0][1]), type="repeat_region")
                feature.qualifiers['standard_name'] = "Repeat_region_5_prime"
                feature.qualifiers['note'] = "poor mapping quality in this region"
                record.features.append(feature)
                RepeatRegionFound_5Prime = True

            if overlap(len(StatDF_Check) - 100, len(StatDF_Check), PoorQualityRegions[-1][0], PoorQualityRegions[-1][1]) >= 1:
                feature = SeqFeature(SimpleLocation(PoorQualityRegions[-1][0], PoorQualityRegions[-1][1]), type="repeat_region")
                feature.qualifiers['standard_name'] = "Repeat_region_3_prime"
                feature.qualifiers['note'] = "poor mapping quality in this region"
                record.features.append(feature)
                RepeatRegionFound_3Prime = True
                    
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
            
            if FeatureName == "putative":
                continue
            
            try:
                FeatureNote = Feature.qualifiers["Note"][0]
            except:
                FeatureNote = Feature.qualifiers["note"][0]
            FeatureSequence = Feature.qualifiers["translation"][0]
            FeatureStart = int(Feature.location.start)
            FeatureEnd = int(Feature.location.end)
            FeatureStrand = Feature.strand
            
            if statistics_file != None:
                               
                NoteAddition = ""
                for PoorQualityRange in PoorQualityRegions:
                    if overlap(FeatureStart,FeatureEnd,PoorQualityRange[0],PoorQualityRange[1]) > 0:
                        FeatureNote += ", poor mapping quality in this region"
            
            feature = SeqFeature(SimpleLocation(FeatureStart, FeatureEnd, FeatureStrand), type="CDS")
            feature.qualifiers['translation'] = FeatureSequence
            feature.qualifiers['note'] = FeatureNote
            feature.qualifiers['product'] = FeatureName
            record.features.append(feature)
            
            with open(ProjectName + "_FeatureTable.txt", "a") as f:
                if FeatureStrand == 1:
                    print(str(FeatureStart+1) + "\t" + str(FeatureEnd) + "\t" + "CDS", file=f)
                if FeatureStrand == -1:
                    print(str(FeatureEnd) + "\t" + str(FeatureStart+1) + "\t" + "CDS", file=f)
                print("\t" + "product " + str(FeatureName), file=f)
                print("\t" + "note " + str(FeatureNote), file=f)
                f.close()

    SeqIO.write(record, ProjectName + "_Annotated.gb", "gb")

    if statistics_file != None:
        def AddRepeatRegionsToFeatureTable(Feature_Table, PoorQualityRegions, RepeatRegionFound_5Prime, RepeatRegionFound_3Prime):
                
            if RepeatRegionFound_5Prime == True:
                Start_5prime = PoorQualityRegions[0][0] + 1
                End_5prime = PoorQualityRegions[0][1]
                with open(Feature_Table,"a") as f:
                    print(str(Start_5prime) + "\t" + str(End_5prime) + "\t" + "repeat_region", file=f)
                    print("\t" + "standard_name " + "Repeat_region_5_prime", file=f)
                    print("\t" + "note " + "poor mapping quality in this region", file=f)
                    f.close()
            
            if RepeatRegionFound_3Prime == True:
                Start_3prime = PoorQualityRegions[-1][0] + 1
                End_3prime = PoorQualityRegions[-1][1]
                with open(Feature_Table,"a") as f:
                    print(str(Start_3prime) + "\t" + str(End_3prime) + "\t" + "repeat_region", file=f)
                    print("\t" + "standard_name " + "Repeat_region_3_prime", file=f)
                    print("\t" + "note " + "poor mapping quality in this region", file=f)
                    f.close()
                

        AddRepeatRegionsToFeatureTable(ProjectName + "_FeatureTable.txt",PoorQualityRegions, RepeatRegionFound_5Prime, RepeatRegionFound_3Prime)
    
    return ProjectName + "_FeatureTable.txt"

def FixGenbankMakeFeatureTable_Legacy(ProjectName, GBKFile, ConcensusFasta, statistics_file):

    def FindPoorQualityRegions(statistics_file, Quality=59, Streak=100):
        
        """
        Find low quality regions of a _statistics file, such as one produced by samtools mpileup. 

        Parameters
        ---
        statistics_file : The path to a .txt _statistics file, or a pd.dataframe constructed as one. Quality in the df.dataframe can be either in Phred Quality scores or their translated versions.
        
        >>> Example
        "21g5i9hj_2024_01_30_statistics.txt"
        
        >>> Example
            Genome          Position  Coverage  Quality
        0   Example_Name    1          4        0.0
        1   Example_Name    2          4        0.0
        2   Example_Name    3          4        0.0
        
        Quality: An integer or float. Maxiumum Average (by base pair) Quality Desired for Region
        
        Streak: A positive integer. Length of positions with low quality in a row before being declared as a low quality.

        Returns
        ---
        A list of two-length lists detailing the boundaries of poorly sequenced regions.

        Example
        ---
        >>> df = FindPoorQualityRegions("21g5i9hj_2024_01_30_statistics.txt")
        """
        if type(statistics_file) is pd.DataFrame:
            try: 
                statistics_file['Quality'] = statistics_file['Quality'].apply(lambda x: round(mean([ord(i) for i in x]), 1) - 33)
            except:
                try:
                    statistics_file['Mapping_Quality'] = statistics_file['Mapping_Quality'].apply(lambda x: round(mean([ord(i) for i in x]), 1) - 33)
                except:
                    pass
        else:
            try:
                statistics_file = pd.read_csv(statistics_file, sep='\t', header = None, names = ['Genome', 'Position', 'Coverage', 'Quality'], quoting=3)
                statistics_file['Quality'] = statistics_file['Quality'].apply(lambda x: round(mean([ord(i) for i in x]), 1) - 33)
            except:
                statistics_file = pd.read_csv(statistics_file, sep='\t', header = None, names = ['Genome', 'Position', 'Coverage', 'Mapping_Quality', 'Quality_Cigar'], quoting=3)
                statistics_file['Mapping_Quality'] = statistics_file['Mapping_Quality'].apply(lambda x: round(mean([ord(i) for i in x]), 1) - 33)

        def getStreaks(L,N,minSize):
            breaks = [i for i,(a,b) in enumerate(zip(L,L[1:]),1) if (a<=N)!=(b<=N)]
            return [[s,e-1] for s,e in zip([0]+breaks,breaks+[len(L)])
                                if e-s>=minSize and L[s]<=N]
        
        try:
            return getStreaks(statistics_file['Quality'], Quality, Streak)
        except:
            return getStreaks(statistics_file['Mapping_Quality'], Quality, Streak)

    PoorQualityRegions = FindPoorQualityRegions(statistics_file)

    #Step1: Add Concensus Data to Record

    for record in SeqIO.parse(ConcensusFasta, "fasta"):
        record.annotations = { "molecule_type" : "DNA" }
        
    with open(ProjectName + "_FeatureTable.txt", "w") as f:
        f.write(">Feature " + record.id + "\n")
        f.close()

    #Step2: Clean up Data Entries
    GBKFile = SeqIO.read(open(GBKFile, "r"), "genbank")
    #Add Repeat Regions to GB
    if statistics_file != None:
        feature = SeqFeature(SimpleLocation(PoorQualityRegions[0][0], PoorQualityRegions[0][1]), type="repeat_region")
        feature.qualifiers['standard_name'] = "Repeat_region_5_prime"
        feature.qualifiers['note'] = "poor mapping quality in this region"
        record.features.append(feature)
        feature = SeqFeature(SimpleLocation(PoorQualityRegions[-1][0], PoorQualityRegions[-1][1]), type="repeat_region")
        feature.qualifiers['standard_name'] = "Repeat_region_3_prime"
        feature.qualifiers['note'] = "poor mapping quality in this region"
        record.features.append(feature)

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
            
            if FeatureName == "putative":
                continue
            
            try:
                FeatureNote = Feature.qualifiers["Note"][0]
            except:
                FeatureNote = Feature.qualifiers["note"][0]
            FeatureSequence = Feature.qualifiers["translation"][0]
            FeatureStart = int(Feature.location.start)
            FeatureEnd = int(Feature.location.end)
            FeatureStrand = Feature.strand
            
            if statistics_file != None:
                
                def overlap(start1, end1, start2, end2):
                    """how much does the range (start1, end1) overlap with (start2, end2)"""
                    return max(max((end2-start1), 0) - max((end2-end1), 0) - max((start2-start1), 0), 0)
                
                NoteAddition = ""
                for PoorQualityRange in PoorQualityRegions:
                    if overlap(FeatureStart,FeatureEnd,PoorQualityRange[0],PoorQualityRange[1]) > 0:
                        FeatureNote += ", poor mapping quality in this region"
            
            feature = SeqFeature(SimpleLocation(FeatureStart, FeatureEnd, FeatureStrand), type="CDS")
            feature.qualifiers['translation'] = FeatureSequence
            feature.qualifiers['note'] = FeatureNote
            feature.qualifiers['product'] = FeatureName
            record.features.append(feature)
            
            with open(ProjectName + "_FeatureTable.txt", "a") as f:
                if FeatureStrand == 1:
                    print(str(FeatureStart+1) + "\t" + str(FeatureEnd) + "\t" + "CDS", file=f)
                if FeatureStrand == -1:
                    print(str(FeatureEnd) + "\t" + str(FeatureStart+1) + "\t" + "CDS", file=f)
                print("\t" + "product " + str(FeatureName), file=f)
                print("\t" + "note " + str(FeatureNote), file=f)
                f.close()

    SeqIO.write(record, ProjectName + "_Annotated.gb", "gb")

    if statistics_file != None:
        def AddRepeatRegionsToFeatureTable(Feature_Table, PoorQualityRegions):
            Start_5prime = PoorQualityRegions[0][0] + 1
            End_5prime = PoorQualityRegions[0][1]
            
            Start_3prime = PoorQualityRegions[-1][0] + 1
            End_3prime = PoorQualityRegions[-1][1]
            
            with open(Feature_Table,"a") as f:
                print(str(Start_5prime) + "\t" + str(End_5prime) + "\t" + "repeat_region", file=f)
                print("\t" + "standard_name " + "Repeat_region_5_prime", file=f)
                print("\t" + "note " + "poor mapping quality in this region", file=f)
                print(str(Start_3prime) + "\t" + str(End_3prime) + "\t" + "repeat_region", file=f)
                print("\t" + "standard_name " + "Repeat_region_3_prime", file=f)
                print("\t" + "note " + "poor mapping quality in this region", file=f)
                f.close()

        AddRepeatRegionsToFeatureTable(ProjectName + "_FeatureTable.txt",PoorQualityRegions)
    
    return ProjectName + "_FeatureTable.txt"

def FindPoorQualityRegions_Legacy(statistics_file, Quality=59, Streak=100):
    
    """
    Find low quality regions of a _statistics file, such as one produced by samtools mpileup. 

    Parameters
    ---
    statistics_file : The path to a .txt _statistics file, or a pd.dataframe constructed as one. Quality in the df.dataframe can be either in Phred Quality scores or their translated versions.
    
    >>> Example
    "21g5i9hj_2024_01_30_statistics.txt"
    
    >>> Example
        Genome          Position  Coverage  Quality
    0   Example_Name    1          4        0.0
    1   Example_Name    2          4        0.0
    2   Example_Name    3          4        0.0
    
    Quality: An integer or float. Maxiumum Average (by base pair) Quality Desired for Region
    
    Streak: A positive integer. Length of positions with low quality in a row before being declared as a low quality.

    Returns
    ---
    A list of two-length lists detailing the boundaries of poorly sequenced regions.

    Example
    ---
    >>> df = FindPoorQualityRegions("21g5i9hj_2024_01_30_statistics.txt")
    """
    if type(statistics_file) is pd.DataFrame:
        try: 
            statistics_file['Quality'] = statistics_file['Quality'].apply(lambda x: round(mean([ord(i) for i in x]), 1) - 33)
        except:
            try:
                statistics_file['Mapping_Quality'] = statistics_file['Mapping_Quality'].apply(lambda x: round(mean([ord(i) for i in x]), 1) - 33)
            except:
                pass
    else:
        try:
            statistics_file = pd.read_csv(statistics_file, sep='\t', header = None, names = ['Genome', 'Position', 'Coverage', 'Quality'], quoting=3)
            statistics_file['Quality'] = statistics_file['Quality'].apply(lambda x: round(mean([ord(i) for i in x]), 1) - 33)
        except:
            statistics_file = pd.read_csv(statistics_file, sep='\t', header = None, names = ['Genome', 'Position', 'Coverage', 'Mapping_Quality', 'Quality_Cigar'], quoting=3)
            statistics_file['Mapping_Quality'] = statistics_file['Mapping_Quality'].apply(lambda x: round(mean([ord(i) for i in x]), 1) - 33)

    def getStreaks(L,N,minSize):
        breaks = [i for i,(a,b) in enumerate(zip(L,L[1:]),1) if (a<=N)!=(b<=N)]
        return [[s,e-1] for s,e in zip([0]+breaks,breaks+[len(L)])
                            if e-s>=minSize and L[s]<=N]
    
    try:
        return getStreaks(statistics_file['Quality'], Quality, Streak)
    except:
        return getStreaks(statistics_file['Mapping_Quality'], Quality, Streak)


#### NEW FUNCTIONS FOR ILLUMINAONLY PIPELINE
def MoveReferenceFile(ReferenceDirectory, PredictedReference, FileType):
    if not os.path.isfile(PredictedReference + FileType):
        shutil.copy2(ReferenceDirectory + PredictedReference + FileType, os.getcwd())
    else:
        pass
    return PredictedReference + FileType

def PredictReference(DeNovo_Scaffold, ProjectName):
    BlastInput =  DeNovo_Scaffold
    BlastReference = "/app/01_References/ReferenceFastas.fa" #docker container file
    BlastOutput = ProjectName + "_predictreference_blastoutput.csv"
    os.system(str(NcbiblastnCommandline(cmd='blastn', query = BlastInput, subject= BlastReference, max_hsps = 1, out = BlastOutput, outfmt = "6 qseqid sseqid evalue pident bitscore length qlen slen qstart qend sstart send", evalue=0.001)))
    BlastDF = pd.read_csv(BlastOutput, sep = '\t', names=('qseqid', 'sseqid', 'evalue', 'pident', 'bitscore', 'length', 'qlen', 'slen', 'qstart', 'qend', 'sstart', 'send'))
    BlastDF = BlastDF[BlastDF['qlen'] >= 1250]
    try:
        PredictedReference = BlastDF.groupby('sseqid')['bitscore'].sum().sort_values(ascending=False).reset_index().iloc[0]['sseqid']
        return PredictedReference
    except:
        return None

#### NEW FUNCTIONS FOR MASTERPIPELINE

def slugify(value, allow_unicode=False):
    """
    Taken from https://github.com/django/django/blob/master/django/utils/text.py
    Convert to ASCII if 'allow_unicode' is False. Convert spaces or repeated
    dashes to single dashes. Remove characters that aren't alphanumerics,
    underscores, or hyphens. Convert to lowercase. Also strip leading and
    trailing whitespace, dashes, and underscores.
    """
    value = str(value)
    if allow_unicode:
        value = unicodedata.normalize('NFKC', value)
    else:
        value = unicodedata.normalize('NFKD', value).encode('ascii', 'ignore').decode('ascii')
    value = re.sub(r'[^\w\s-]', '', value.lower())
    return re.sub(r'[-\s]+', '-', value).strip('-_')

def MiniMapScaffoldToReference(GenomeFile, ScaffoldFile, ProjectName, Suffix):
    ScaffoldMapped_SAM = ProjectName + "_" + Suffix + ".sam"
    os.system("minimap2 -ax asm5 " + GenomeFile + " " + ScaffoldFile + " > " + ScaffoldMapped_SAM)
    print("***_MiniMap  Complete_***")
    return ScaffoldMapped_SAM

def Empty_File_Check(file):
    if file != None:
        if os.path.exists(file):
            File_Size = os.stat(file).st_size
            if(File_Size == 0):
                print(file + " has zero data  (Empty_File_Check)")
                return True
            else:
                return False
        else:
            print(file + " does not exist (Empty_File_Check)")
            return True
    else:
        print(file + " is an empty variable  (Empty_File_Check)")
        return True

def Contig_Extender_for_Polish_V02(ProjectName, MinLengthThreshold, MaxLengthThreshold, MaxSubSample, Threshold, Require_Multiple, MinionMappedFirst10, MinionMappedLast10, First_Output_Fasta, Last_Output_Fasta, First_Output_ALN, Last_Output_ALN):
    # FileOutputs
    x5_Prime_NonPolish_FlankingFile = ProjectName + "_5_prime_non_polish.fa"
    x3_Prime_NonPolish_FlankingFile = ProjectName + "_3_prime_non_polish.fa"
    # Import SAM files as Pandas, first 11 colulmns are mandatory, drop unessary columns
    colnames = ['Qname','Flag','Rname','Pos','MapQ','Cigar','RNext','PNExt','Tlen','Seq','Qual']
    droppedColumns = ['Rname','MapQ','RNext','PNExt','Tlen','Qual']
    if Empty_File_Check(MinionMappedFirst10) == False:
        First_RAW_Table = pd.read_csv(MinionMappedFirst10, sep="\t", quoting=3, usecols=range(11), header = None, index_col=False)
        Cancel_5_extend = False
    else:
        Cancel_5_extend = True
        First_RAW_Table = pd.DataFrame(columns = colnames)
        print("5' region will not be extended")
    if Empty_File_Check(MinionMappedLast10) == False:
        Last_RAW_Table = pd.read_csv(MinionMappedLast10, sep="\t", quoting=3, usecols=range(11), header = None, index_col=False)
        Cancel_3_extend = False
    else:
        Cancel_3_extend = True
        Last_RAW_Table = pd.DataFrame(columns = colnames)
        print("3' region will not be extended")
    
    First_RAW_Table.columns = colnames
    Last_RAW_Table.columns = colnames
    #
    MinReadNumberForAlignment = 12
    # DropColumns that are not needed
    First_Working_Table = First_RAW_Table.drop(columns=droppedColumns)
    Last_Working_Table = Last_RAW_Table.drop(columns=droppedColumns)
    #
    #########################################################################################################
    #                               Explanation of Regex                                                    #
    # Filter by Regex:                                                                                      #
    #       Regrex for 5' Region of Genome = '^[0-9]{3,4}S'                                                 #
    #              "^" = StartofString "[0-9]" = AnyNumber "{3,4}" Appearing 3 to 4 times "S"               #
    #       Regrex for 3' Region of Genome = '[0-9]{3,4}S$'                                                 #
    #              "[0-9]" = AnyNumber "{3,4}" Appearing 3 to 4 times "S" $= EndofString                    #
    #                                                                                                       #
    #########################################################################################################
    #
    First_Working_Table = First_Working_Table[First_Working_Table.Cigar.str.contains('^[0-9]{3,4}S')].reindex()
    Last_Working_Table = Last_Working_Table[Last_Working_Table.Cigar.str.contains('[0-9]{3,4}S$')].reindex()
    #
    # Determine Flanking Nucleotides by Spliting @ 'S'
    if Cancel_5_extend != True:
        SplitSampleName_First = First_Working_Table.Cigar.str.split(pat = 'S', expand = True)[0]
    #
    ########      Overly complicated method for 3' end because [-1] would not work      #####################
    if Cancel_3_extend != True:
        # write the string mirrored                                                                         #
        Last_Working_Table['Cigar'] = Last_Working_Table['Cigar'].apply(lambda x: x[::-1])                  #
        # replace all letters besides S with #                                                              #
        Last_Working_Table['Cigar'] = Last_Working_Table.Cigar.str.replace('[A-RT-Z]', '#', regex = True)   #
        # Split @ First #                                                                                   #
        Last_Working_Table['Cigar'] = Last_Working_Table.Cigar.str.split(pat = '#', expand = True)[0]       #
        # Correct Mirror                                                                                    #
        Last_Working_Table['Cigar'] = Last_Working_Table['Cigar'].apply(lambda x: x[::-1])                  #
        # Determine Flanking Nucleotides by Spliting @ 'S'                                                  #
        SplitSampleName_Last = Last_Working_Table.Cigar.str.split(pat = 'S', expand = True)[0]              #
    #                                                                                                       #
    #########################################################################################################
    #
    # Setting Up Bin and Count the length of the flanking sequences Variables
    bins = range(MinLengthThreshold, MaxLengthThreshold, 500)
    labels = bins[:-1]

    # 5' end of genome
    if Cancel_5_extend != True:
        # Create New Column and Convert to Numeric
        First_Working_Table['Flanking'] = SplitSampleName_First
        First_Working_Table['Flanking'] = pd.to_numeric(First_Working_Table['Flanking'])
        # Calculate the length of Sequence
        First_Working_Table['len'] = First_Working_Table['Seq'].str.len()
        # Extract First N Nucleotides(flanking 5' end)
        First_Working_Table['Seq2'] = First_Working_Table.apply(lambda x: x['Seq'][: x['Flanking'] + 1], 1)
        # Calculate the length of new Sequence
        First_Working_Table['len2'] = First_Working_Table['Seq2'].str.len()
        # Filter Reads above certain thereshold
        First_Working_Table = First_Working_Table[First_Working_Table.len2 >= MinLengthThreshold].reindex()
        # Bin
        First_Bin_Count = First_Working_Table.groupby(pd.cut(First_Working_Table['len2'], bins=bins, labels=labels)).size().to_frame()
        First_Bin_Count.insert(0, 'LowerBoundary', First_Bin_Count.index)
        First_Bin_Count = First_Bin_Count.rename(columns={0: "Count"})
        # only include bins that contain > MinReadNumberForAlignment
        First_Bin_Count = First_Bin_Count[First_Bin_Count.Count >= MinReadNumberForAlignment].reset_index(drop=True)
        
        # only perform this process if there are enough bins to justify: min length = 1000-1500
        if len(First_Bin_Count) > 1:
            First_Smallest_Size = First_Bin_Count.LowerBoundary[len(First_Bin_Count)- 1]
            First_Largest__Size = First_Smallest_Size + 500
            First_Working_Table = First_Working_Table[First_Working_Table.len2 >= First_Smallest_Size].reindex()
            First_Working_Table = First_Working_Table[First_Working_Table.len2 <= First_Largest__Size].reindex()
            Cancel_5_extend = False
        else:
            Cancel_5_extend = True
            print("5' region will not be extended")
    
    # 3' end of genome
    if Cancel_3_extend != True:
        # Create New Column and Convert to Numeric
        Last_Working_Table['Flanking'] = SplitSampleName_Last
        Last_Working_Table['Flanking'] = pd.to_numeric(Last_Working_Table['Flanking'])
        # Calculate the length of Sequence
        Last_Working_Table['len'] = Last_Working_Table['Seq'].str.len()
        # Extract  Last N nucleotides (flanking 3' end)
        Last_Working_Table['Seq2'] = Last_Working_Table.apply(lambda x: x['Seq'][x['len'] - x['Flanking'] : x['len'] + 1], 1)
        # Calculate the length of new Sequence
        Last_Working_Table['len2'] = Last_Working_Table['Seq2'].str.len()
        # Filter Reads above certain thereshold
        Last_Working_Table = Last_Working_Table[Last_Working_Table.len2 >= MinLengthThreshold].reindex()
        # Bin
        Last__Bin_Count = Last_Working_Table.groupby(pd.cut(Last_Working_Table['len2'], bins=bins, labels=labels)).size().to_frame()
        Last__Bin_Count.insert(0, 'LowerBoundary', Last__Bin_Count.index)
        Last__Bin_Count = Last__Bin_Count.rename(columns={0: "Count"})
        # only include bins that contain > MinReadNumberForAlignment
        Last__Bin_Count = Last__Bin_Count[Last__Bin_Count.Count >= MinReadNumberForAlignment].reset_index(drop=True)
        if len(Last__Bin_Count) > 1:
            Last__Smallest_Size = Last__Bin_Count.LowerBoundary[len(Last__Bin_Count)- 1]
            Last__Largest__Size = Last__Smallest_Size + 500
            Last_Working_Table = Last_Working_Table[Last_Working_Table.len2 >= Last__Smallest_Size].reindex()
            Last_Working_Table = Last_Working_Table[Last_Working_Table.len2 <= Last__Largest__Size].reindex()
            Cancel_3_extend = False
        else:
            Cancel_3_extend = True
            print("3' region will not be extended")
    #
    #
    # Extract only a subset of reads to save Computational Time (if needed)
    First_Working_Table_row_Count = First_Working_Table.shape[0]
    Last_Working_Table_row_Count = Last_Working_Table.shape[0]
    if First_Working_Table_row_Count > MaxSubSample:
        First_Working_Table = First_Working_Table.sample(n = MaxSubSample, random_state=1)
    if Last_Working_Table_row_Count > MaxSubSample:
        Last_Working_Table = Last_Working_Table.sample(n = MaxSubSample, random_state=1)
    #
    # Print to CSV (testing purposes only)
    #First_Working_Table.to_csv(First_Output_Fasta[:-2]+"csv", sep=",")
    #Last_Working_Table.to_csv(Last_Output_Fasta[:-2]+"csv", sep=",")
    #
    # Print Flanking Reads to Fasta Files
    if Cancel_5_extend == False:
        First_Name_List = First_Working_Table['Qname'].tolist()
        First_Seq_List = First_Working_Table['Seq2'].tolist()
        with open(First_Output_Fasta, "w") as f:
            for Name, Sequence in zip(First_Name_List, First_Seq_List):
                f.write(">" + Name + "\n" + Sequence + "\n")
        f.close()  
    #
    if Cancel_3_extend == False:
        Last_Name_List = Last_Working_Table['Qname'].tolist()
        Last_Seq_List = Last_Working_Table['Seq2'].tolist()
        with open(Last_Output_Fasta, "w") as f:
            for Name, Sequence in zip(Last_Name_List, Last_Seq_List):
                f.write(">" + Name + "\n" + Sequence + "\n")
        f.close()
    
    #
    print("***_Contig Extender 1/6 Complete_***")
    #
    # Create Alignment with muscle
    if Cancel_5_extend == False:
        os.system("muscle -in " + First_Output_Fasta + " -out " + First_Output_ALN)
    print("***_Contig Extender 2/6_Complete_***")
    if Cancel_3_extend == False:
        os.system("muscle -in " + Last_Output_Fasta + " -out " + Last_Output_ALN)
    print("***_Contig Extender 3/6 Complete_***")
    #
    # Extract Simple Consensus from alignments
    if Cancel_5_extend == False:
        alignment_5prime  = AlignIO.read(open(First_Output_ALN), "fasta")
        consensus_5prime = AlignInfo.SummaryInfo(alignment_5prime).gap_consensus(threshold=Threshold, ambiguous='N', require_multiple=Require_Multiple)
        consensus_5prime = str(consensus_5prime).replace('-','')
    else:
        consensus_5prime = "NNNNN"
    print("***_Contig Extender 4/6 Complete_***")
    #
    if Cancel_3_extend == False:
        alignment_3prime  = AlignIO.read(open(Last_Output_ALN), "fasta")
        consensus_3prime = AlignInfo.SummaryInfo(alignment_3prime).gap_consensus(threshold=Threshold, ambiguous='N', require_multiple=Require_Multiple)
        consensus_3prime = str(consensus_3prime).replace('-','')
    else:
        consensus_3prime = "NNNNN"
    print("***_Contig Extender 5/6 Complete_***")
    #
    # Create Final Fasta
    with open(x5_Prime_NonPolish_FlankingFile, "w") as f:
        f.write(">" + ProjectName + "5_Prime_Flank_nonpolished\n" + str(consensus_5prime) + "\n")
    f.close()
    with open(x3_Prime_NonPolish_FlankingFile, "w") as f:
        f.write(">" + ProjectName + "3_Prime_Flank_nonpolished\n" + str(consensus_3prime) + "\n")
    f.close()
    print("***Contig Extender 6/6_Complete_***")
    return x5_Prime_NonPolish_FlankingFile, x3_Prime_NonPolish_FlankingFile

def SpadesDeNovo(ExtractedRead01, ExtractedRead02, ProjectName, CombinedMinionReadFile = None, ExtractedRead03 = None, ExtractedRead04 = None, threads = "3", RAM = "10", Kmers = "21,33,55,75,99"):
    OutputDirectory = ProjectName + "_DeNovoAssembly"
    Assembled_Scaffold_File = OutputDirectory + "/scaffolds.fasta"
    
    if CombinedMinionReadFile == None:
        MinionString = ""
    else: 
        MinionString = " --nanopore " + CombinedMinionReadFile
        print("***Running De Novo Assembly with Nanopore Reads:" + CombinedMinionReadFile + "***")
    
    if ExtractedRead03 == None or ExtractedRead04 == None:
        Read_String = "-1 " + ExtractedRead01 + " -2 " + ExtractedRead02
        print("***Running De Novo Assembly with Illumina Reads:" + ExtractedRead01 + ExtractedRead02 + "***")
    else:
        Read_String = "--pe1-1 " + ExtractedRead01 + " --pe1-2 " + ExtractedRead02 + " --pe2-1 " + ExtractedRead03 + " --pe2-2 " + ExtractedRead04
        print("***Running De Novo Assembly with Illumina Reads:" + ExtractedRead01 + ExtractedRead02 + ExtractedRead03 + ExtractedRead04 + "***")
    
    #os.system("spades.py --careful " + Read_String + MinionString + " -t " + str(threads) + " -m " + str(RAM) + " -k " + str(Kmers) + " -o " + OutputDirectory)
    os.system("spades.py --careful " + Read_String + MinionString + " -t " + str(threads) + " -k " + str(Kmers) + " -o " + OutputDirectory)
    
    
    print("***De Novo Assembly_Complete_***")
    
    return Assembled_Scaffold_File

def RemoveSwineContigs(DeNovo_Scaffold, ProjectName):
    BlastInput =  DeNovo_Scaffold
    BlastReference = "/app/06_Sus_scrofa_Reference/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa" #docker container file
    BlastOutput = ProjectName + "_predict_swine_contigs_blastoutput.csv"
    ASFV_only_Scaffold = ProjectName + "_scaffolds.fasta"
    #
    print("***Checking for Swine Contigs in Assembly***")
    os.system(str(NcbiblastnCommandline(cmd='blastn', query = BlastInput, subject= BlastReference, max_hsps = 1, max_target_seqs = 1, out = BlastOutput, outfmt = "6 qseqid sseqid evalue pident bitscore length qlen slen qstart qend sstart send", evalue=0.001)))
    #
    if Empty_File_Check(BlastOutput) == True:
        # Move and Rename Scaffold file for retention
        print('No Swine Contigs')
        MoveReferenceFile(ReferenceDirectory="" + ProjectName + "_DeNovoAssembly/", PredictedReference = "scaffolds", FileType=".fasta")
        Rename_File(OldFile = "scaffolds.fasta", NewFile = ASFV_only_Scaffold)
        return ASFV_only_Scaffold
    else:
        print('***Removing Swine Contigs***')
        
        BlastDF = pd.read_csv(BlastOutput, sep = '\t', names=('qseqid', 'sseqid', 'evalue', 'pident', 'bitscore', 'length', 'qlen', 'slen', 'qstart', 'qend', 'sstart', 'send'), usecols=['qseqid'])
        Contig_List = BlastDF['qseqid'].tolist()
        with open(ASFV_only_Scaffold, "w") as f:
            for seq_record in SeqIO.parse(DeNovo_Scaffold, "fasta"):
                if not seq_record.id in Contig_List:
                    f.write(">" + seq_record.id  + "\n" + str(seq_record.seq) + "\n")
        f.close()
        return ASFV_only_Scaffold


def MiniMapToReferenceSimple2(GenomeFileTemp, CombinedMinionReadFile, ProjectName, Suffix):
    MinionMappedTempSAM = ProjectName + "_" + Suffix + "_Minion_Mapped.sam"
    MinionMappedTempBAM = ProjectName + "_" + Suffix + "_Minion_Mapped.bam"
    MinionMappedTempSOR = ProjectName + "_" + Suffix + "___Minion_Sort.sam"
    os.system("minimap2 -ax map-ont " + GenomeFileTemp + " " + CombinedMinionReadFile + " > " + MinionMappedTempSAM)
    print("***_MiniMap 1/3 Complete_***")
    os.system("samtools view -b -F 4 " + MinionMappedTempSAM + " > " + MinionMappedTempBAM)
    print("***_MiniMap 2/3 Complete_***")
    os.system("samtools sort " + MinionMappedTempBAM + " -o " + MinionMappedTempSOR)
    print("***_MiniMap 3/3 Complete_***")
    return [MinionMappedTempSAM, MinionMappedTempBAM, MinionMappedTempSOR]


#### Legacy Functions
'''
def MapToReference_Tweaked_K(Trimmed_Read01, Trimmed_Read02, Lane, Threads, ReferenceDirectory, Reference, ProjectName, Seed):
    if Trimmed_Read01 == None or Trimmed_Read02 == None:
        OutputLocation = ""
        return OutputLocation
    
    ReferenceLocation = ReferenceDirectory + Reference
    OutputLocation = ProjectName + "_mapped_" + Reference.replace(".fa","") + "_Lane" + str(Lane) + ".bam"
    os.system("bwa-mem2 mem -t " + Threads + " -k " + Seed + " " + ReferenceLocation + " " + Trimmed_Read01 + " " + Trimmed_Read02 + " | samtools sort -o " + OutputLocation + " " + "-")
    
    print("***_Illumina Mapping Complete_***")
    return OutputLocation

def SpadesDeNovo_Legacy_02(ExtractedRead01, ExtractedRead02, ProjectName, CombinedMinionReadFile = None, threads = "3", RAM = "10", Kmers = "21,33,55,75,99"):
    #RAM is in Gigabytes
    OutputDirectory = ProjectName + "_DeNovoAssembly"
    Assembled_Scaffold_File = OutputDirectory + "/scaffolds.fasta"
    if CombinedMinionReadFile == None:
        os.system("spades.py --careful -1 " + ExtractedRead01 + " -2 " + ExtractedRead02 + " -t " + str(threads) + " -m " + str(RAM) + " -k " + str(Kmers) + " -o " + OutputDirectory)
        print("***De Novo Assembly_Complete_***")
        return Assembled_Scaffold_File
    else:
        os.system("spades.py --careful -1 " + ExtractedRead01 + " -2 " + ExtractedRead02 + " --nanopore " + CombinedMinionReadFile + " -t " + str(threads) + " -m " + str(RAM) + " -k " + str(Kmers) + " -o " + OutputDirectory)
        print("***De Novo Assembly_Complete_***")
        return Assembled_Scaffold_File

def SpadesDeNovo_Legacy(ExtractedRead01, ExtractedRead02, CombinedMinionReadFile, ProjectName):
    os.system("spades.py --careful -1 " + ExtractedRead01 + " -2 " + ExtractedRead02 + " --nanopore " + CombinedMinionReadFile + " -t 3 -m 10 -k 21,33,55,75,99 -o " + ProjectName)
    Assembled_Scaffold_File = ProjectName + "/scaffolds.fasta"
    print("***De Novo Assembly_Complete_***")
    return Assembled_Scaffold_File

def Move_File_Legacy(Date, ProjectName, List_of_Files):
    NewDirectory = Date + "_Project_" + ProjectName +"/"
    if not os.path.exists(NewDirectory):
        os.makedirs(NewDirectory)
    else:
        pass
    for file in List_of_Files:
        if file != None:
            if os.path.exists(file):
                shutil.move(file, NewDirectory + file)
        else:
            print("empty variable")

def Make_Directory_Legacy(ProjectName):
    NewDirectory = ProjectName + "_BLASTOUTPUTBANK/"
    if not os.path.exists(NewDirectory):
        os.makedirs(NewDirectory)
        return NewDirectory
    else:
        return NewDirectory

def ExtractFromMerged_UnMappedReads_TEST(Input, Reference, ProjectName):
    
    # Explanation:
    #   -b  =   Output in the BAM format. 
    #   -f  =   Only output alignments with all bits set in INT present in the FLAG field
    #   -F  =   Do not output alignments with any bits set in INT present in the FLAG field
    #   0x2 = 	PROPER_PAIR	each segment properly aligned according to the aligner
    #   0x4 =	UNMAP, segment unmapped
    os.system("samtools view -b -f12 " + Input + " | samtools sort -n -o " + ProjectName + "_mappedReads_" + Reference + "_sort.bam")
    print("***_Reads Extracted Complete 1/2_***")

    #extract as fastq-paired Reads
    ExtractedRead01 = ProjectName + "_" + Reference + "_TrimmedExtracted1.fastq.gz"
    ExtractedRead02 = ProjectName + "_" + Reference + "_TrimmedExtracted2.fastq.gz"
    
    os.system("samtools fastq -@ 8 " + ProjectName + "_mappedReads_" + Reference + "_sort.bam -1 " + ExtractedRead01 + " -2 " + ExtractedRead02 + " -0 /dev/null -s /dev/null -n")
    
    print("***_Reads Extracted Complete 2/2_***")
    return ExtractedRead01, ExtractedRead02
'''
