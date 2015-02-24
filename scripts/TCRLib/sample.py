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
	self.TCRcaller = TCRcaller

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
	for filePairId,readCount,fastq1,fastq2,sampleId in self.TCRcaller.database.getFastqs():
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

	import glob
	import os
	import sys
	import subprocess
	import multiprocessing

        # find script path
        path = os.path.abspath(sys.argv[0])
        path = os.path.abspath('/'.join(sys.argv[0].split('/')[:-2]))

	latestref = newest = max(glob.iglob(path+'/references/ncbi.*.fa'), key=os.path.getctime)

	for fastq in self.getFastqs():
	    
	    sys.stderr.write('aligning '+self.name+' fastq '+fastq[2]+'...')
	    command = [path+'/bin/bowtie2-2.2.4/'+'bowtie2','-x',latestref,'-U',fastq[2],'-U',fastq[3],'-S',self.dataPath+'/'+str(fastq[0])+'.sam','--very-sensitive-local','-p',str(multiprocessing.cpu_count())]
	    #sys.stderr.write('running: '+' '.join(command)+'\n')
	    bt2 = subprocess.Popen(command,stdout=subprocess.PIPE,stderr=subprocess.PIPE )
	    errdata = bt2.communicate()
	    if bt2.returncode != 0:
	        print 'cmd: '+' '.join( command )
	        print 'bt2 view Error code', bt2.returncode, errdata
	        sys.exit()
	    sys.stderr.write('done\n')	    

	    sys.stderr.write('converting '+self.name+' samfile '+str(fastq[0])+' to bam...')
	    bam = open(self.dataPath+'/'+str(fastq[0])+'.bam','wb')
	    command = [path+'/bin/samtools-1.2/'+'samtools','view','-Sb',self.dataPath+'/'+str(fastq[0])+'.sam']
	    samtools = subprocess.Popen(command,stdout=bam,stderr=subprocess.PIPE )
	    errdata = samtools.communicate()
	    if samtools.returncode != 0:
	        print 'cmd: '+' '.join( command )
	        print 'bt2 view Error code', samtools.returncode, errdata
	        sys.exit()
	    sys.stderr.write('done\n')
	    os.remove(self.dataPath+'/'+str(fastq[0])+'.sam')
	    bam.close()

	    sys.stderr.write('sorting '+self.name+' bamfile '+str(fastq[0])+'...')
	    command = [path+'/bin/samtools-1.2/'+'samtools','sort',self.dataPath+'/'+str(fastq[0])+'.bam',self.dataPath+'/'+str(fastq[0])+'.sorted']
	    samtools = subprocess.Popen(command,stdout=subprocess.PIPE,stderr=subprocess.PIPE )
	    errdata = samtools.communicate()
	    if samtools.returncode != 0:
	        print 'cmd: '+' '.join( command )
	        print 'bt2 view Error code', samtools.returncode, errdata
	        sys.exit()
	    sys.stderr.write('done\n')
	    os.remove(self.dataPath+'/'+str(fastq[0])+'.bam')

	    sys.stderr.write('indexing '+self.name+' bamfile '+str(fastq[0])+'...')
	    command = [path+'/bin/samtools-1.2/'+'samtools','index',self.dataPath+'/'+str(fastq[0])+'.sorted.bam']
	    samtools = subprocess.Popen(command,stdout=subprocess.PIPE,stderr=subprocess.PIPE )
	    errdata = samtools.communicate()
	    if samtools.returncode != 0:
	        print 'cmd: '+' '.join( command )
	        print 'bt2 view Error code', samtools.returncode, errdata
	        sys.exit()
	    sys.stderr.write('done\n')

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