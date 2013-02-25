import os
import fnmatch
import subprocess
import sys
import re


def bams_in_one_dir(bams):
    print(bams)
    
    for bam in bams:
        prefix = bam.split("-")[0]
        
        outputdir =   os.curdir + '/'+prefix+"_DepthOfCoverageWithGeneName"   
        os.mkdir(outputdir)
        JavaParm = "java -Xmx4g -jar /miseqdata/tools/GATK/GenomeAnalysisTK.jar "
        GATKAnalysis = "-T DepthOfCoverage "
        InputFile = "-I "+ os.curdir +"/" + bam
        RefFile = " -R /miseqdata/safer/genomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/ucsc.hg19.fasta" 
        Output = " -o "+outputdir+'/'+ prefix
        TargetRegion = " -L " + sys.argv[1] 
        Options = " -omitBaseOutput --summaryCoverageThreshold 30"
    #GeneCoverage = "--calculateCoverageOverGenes /miseqdata/"
        CompleteJavaString = JavaParm+GATKAnalysis+InputFile+RefFile +Output+TargetRegion +Options#+GeneCoverage
        proc = subprocess.Popen(CompleteJavaString,shell=True);
        proc.wait()
        PrintReport(Output)

def bams_in_seperate_dirs():
    dirs = [i for i in os.listdir(os.curdir) if os.path.isdir(i)]

    for dirr in dirs:
        print (os.curdir+'/'+dirr)
        bamfile = [file for file in os.listdir(os.curdir+'/'+dirr) if fnmatch.fnmatch(file,'*recalibrated.bam')]
        print(bamfile)
    
    
    prefix = bamfile[0].split(".")[0]
    print(prefix)
    outputdir =   os.curdir + '/'+dirr+ "_DepthOfCoverageWithGeneName"   
    #os.mkdir(outputdir)
    JavaParm = "java -Xmx4g -jar /miseqdata/tools/GATK/GenomeAnalysisTK.jar "
    GATKAnalysis = "-T DepthOfCoverage "
    InputFile = "-I "+ os.curdir +"/" + bamfile[0]
    RefFile = " -R /miseqdata/safer/genomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/ucsc.hg19.fasta" 
    Output = " -o "+outputdir+'/'+ prefix
    TargetRegion = " -L "+ sys.argv[1]
    Options = " -omitBaseOutput --summaryCoverageThreshold 30"
    GeneCoverage = "--calculateCoverageOverGenes /miseqdata/"
    CompleteJavaString = JavaParm+GATKAnalysis+InputFile+RefFile +Output+TargetRegion +Options#+GeneCoverage
    proc = subprocess.Popen(CompleteJavaString,shell=True);
    proc.wait()
    PrintReport(Output)
#os.chdir("/miseqdata/121110_M00386_0007_000000000-A26PP/Data/Intensities/BaseCalls/Alignment")

def PrintReport(IntervalSummaryFile):
    
    AnnotationDict = CreateAnnotationDict()
    IntervalSummaryFile = re.sub("-o","",IntervalSummaryFile)
    IntervalSummaryFile += '.sample_interval_summary'
    #print (AnnotationDict.get('chrX:153363060-153363188',"FALSE"))
    OutputFile = IntervalSummaryFile.rsplit("/",1)[:-1]
    OutputFileStr = OutputFile[0].strip()+"/DepthOfCoverageAnalysis.txt"
    #OutputFileStr = OutputFileStr.strip()
    OUTFH = open(OutputFileStr,'w')
    with open(IntervalSummaryFile.strip(),'r') as ISF:
        for line in ISF:
            if re.match("Target",line):
                line = re.sub(r"\n","",line)
                OUTFH.write(line+"\tGene Name and Exon\n")
            else:
                line = re.sub(r"\n","",line)
                OUTFH.write("{line}\t{GeneName}\n".format(line=line,GeneName=
                AnnotationDict.get(line.split()[0],"FALSE")))
    
def CreateAnnotationDict():
    AnnotationFile ='/miseqdata/Annotation_files/52GenesWithExonCoordinatesUSEFORDepthOfCovarage.txt'
    AnnotationDict = {}
    with open(AnnotationFile,"r") as AnnFile:
        for line in AnnFile:
            LineList = line.split("\t")
            #Must add one because value are bed format but GATK is 1-based
            AnnotationDict[LineList[1]+":"+str((int(LineList[2])+1))+"-"+LineList[3]] = LineList[0]+" Exon "+LineList[4]
    return AnnotationDict 
    
os.chdir('/miseqdata/130216_M00386_0017_000000000-A35M5/Data/Intensities/BaseCalls/Alignment/')
bams= [bam for bam in os.listdir(os.curdir) if os.path.isfile(os.curdir+"/"+bam) and fnmatch.fnmatch(bam,"*.bam")]

if bams != None:
    bams_in_one_dir(bams)
else:
    bams_in_seperate_dirs()
    
