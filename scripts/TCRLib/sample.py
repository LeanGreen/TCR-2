#! /bin/env python

def main(): pass

class Sample(object):

    def __init__(self, TCRcaller, sampleName=None,sampleId=None,refType=None):
	self.name = sampleName
	self.id = int(sampleId)
	self.refType = refType
	self.path = TCRcaller.path+'/samples/'+self.name
	self.scriptPath = TCRcaller.path+'/samples/'+self.name+'/script'
	self.dataPath   = TCRcaller.path+'/samples/'+self.name+'/data'
	self.logPath    = TCRcaller.path+'/samples/'+self.name+'/logs'
	self.fastqcPath = TCRcaller.path+'/samples/'+self.name+'/fastQC'
	self.dependencies = {}

    @property
    def readCount(self, ):
	tmpCounter = 0
	for filePairId,readCount,fastq1,fastq2,sampleId in AnalysisPipe.database.getFastqs():
	    if int(sampleId) == self.id:
		try: tmpCounter+= readCount
		except TypeError: return 'Unknown'
	return tmpCounter

    def getFastqs(self,):
	self.fastqIds = []
	for filePairId,readCount,fastq1,fastq2,sampleId in AnalysisPipe.database.getFastqs():
	    if int(sampleId) == self.id:
		self.fastqIds.append(filePairId)
		yield [filePairId,readCount,fastq1,fastq2,sampleId]

    def createDirs(self):
        import os
        try: os.makedirs(self.path)
        except OSError:pass
        try: os.makedirs(self.scriptPath)
        except OSError:pass
        try: os.makedirs(self.dataPath)
        except OSError:pass
        try: os.makedirs(self.fastqcPath)
        except OSError:pass
        try: os.makedirs(self.logPath)
        except OSError:pass

    def align(self,):
	#bowtie2 -x references/ncbi.T_cell_receptor_beta.Human.mRNA.fa -U data/mnt/davidson/sandberglab/pipeline2.0/rnaseq/hsa/simone_picelli_T-cells/rawdata/$i/Run00221_L7_1_140213_SN893_0221_BC3H8JACXX_*.fastq.gz -U data/mnt/davidson/sandberglab/pipeline2.0/rnaseq/hsa/simone_picelli_T-cells/rawdata/$i/Run00223_L3_1_140220_SN893_0223_BC3JDHACXX_*.fastq.gz -S $i.sam --very-sensitive-local -p16
	#samtools view -Sb $i.sam > $i.bam
	#samtools sort $i.bam $i.sorted
        #samtools index $i.sorted.bam
	pass
    
    def getAlignedReadsPairs(self, ):
	#samtools view -F 4 -q 20 $i.sorted.bam | cut -f 1,10,11 | awk '{print "@" $1 "\n" $2 "\n+\n" $3}' > $i.mapqMoreThan20.fq
	pass

    def assemble(self,):
#	velveth $i 21 -fastq -short $i.mapqMoreThan20.fq
#	velvetg $i/ -cov_cutoff auto
	pass

    def report(self,):
#	cat $i/contigs.fa
#	mv -v $i* $i
#	echo -e "\n### The reads mapped to the following gi: "
#	samtools view -F 4 -q 20 $i/$i.sorted.bam | cut -f3 | sort | uniq -c;
#	echo -e "\n### same as to the following sequences: "
#	echo -e $(for name in $(samtools view -F 4 -q 20 $i/$i.sorted.bam | cut -f3 | sort | uniq); do grep $name references/ncbi.T_cell_receptor_beta.Human.mRNA.fa;done) | sed s/\>/\\\n/g | grep -P "T|cell|receptor|beta|$" --color=Always
	pass
    

if __name__ == "__main__": main()