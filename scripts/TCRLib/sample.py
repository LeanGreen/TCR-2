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

	latestref = newest = max(glob.iglob(path+'/references/*.fa'), key=os.path.getctime)

	for fastq in self.getFastqs():
	    
	    sys.stderr.write('aligning '+self.name+' fastq '+fastq[2]+'...')
	    if fastq[3]: command = [path+'/bin/bowtie2-2.2.4/'+'bowtie2','-x',latestref,'-U',fastq[2],'-U',fastq[3],'-S',self.dataPath+'/'+str(fastq[0])+'.sam','--very-sensitive-local','-p',str(multiprocessing.cpu_count())]
	    else:        command = [path+'/bin/bowtie2-2.2.4/'+'bowtie2','-x',latestref,'-U',fastq[2],'-S',self.dataPath+'/'+str(fastq[0])+'.sam','--very-sensitive-local','-p',str(multiprocessing.cpu_count())]
	    #sys.stderr.write('running: '+' '.join(command)+'\n')
	    bt2 = subprocess.Popen(command,stdout=subprocess.PIPE,stderr=subprocess.PIPE )
	    errdata = bt2.communicate()
	    if bt2.returncode != 0:
	        print 'cmd: '+' '.join( command )
	        print 'bt2 view Error code', bt2.returncode, errdata
	        sys.exit()
	    sys.stderr.write('done\n')	    

    def postAlign(self,):

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
    
    def getAlignedReadsPairs(self, ):

	import glob
	import os
	import sys
	import subprocess
	import multiprocessing

        # find script path
        path = os.path.abspath(sys.argv[0])
        path = os.path.abspath('/'.join(sys.argv[0].split('/')[:-2]))

	latestref = newest = max(glob.iglob(path+'/references/ncbi.*.fa'), key=os.path.getctime)

	tmpreads = open(self.dataPath+'/mappedReads.tmp','w')
	for fastq in self.getFastqs():
	    sys.stderr.write('getting aligned reads '+self.name+' bamfile '+str(fastq[0])+'...')
	    command = [path+'/bin/samtools-1.2/'+'samtools','view','-F','4','-q','20',self.dataPath+'/'+str(fastq[0])+'.sorted.bam']
	    samtools = subprocess.Popen(command,stdout=tmpreads,stderr=subprocess.PIPE )
	    errdata = samtools.communicate()
	    if samtools.returncode != 0:
	        print 'cmd: '+' '.join( command )
	        print 'samtools view Error code', samtools.returncode, errdata
	        sys.exit()
	    sys.stderr.write('done\n')
	    os.remove(self.dataPath+'/'+str(fastq[0])+'.sorted.bam')
	    os.remove(self.dataPath+'/'+str(fastq[0])+'.sorted.bam.bai')
	tmpreads.close()

	tmpreads2 = open(self.dataPath+'/mappedReads.2.tmp','w')
	sys.stderr.write('converting aligned reads to fq for '+self.name+'...')
	command = ['cut','-f','1,10,11',tmpreads.name]
	cut = subprocess.Popen(command,stdout=tmpreads2,stderr=subprocess.PIPE )
	errdata = cut.communicate()
	if cut.returncode != 0:
	    print 'cmd: '+' '.join( command )
	    print 'cut view Error code', cut.returncode, errdata
	    sys.exit()
	tmpreads2.close()
	os.remove(tmpreads.name)

	output = open(self.dataPath+'/mapqMoreThan20.fq','w')
	with open(tmpreads2.name) as infile:
	    for line in infile:
		line=line.split()
		output.write('@'+line[0]+'\n'+line[1]+'\n+\n'+line[2]+'\n')
	sys.stderr.write('done\n')
	output.close()
	os.remove(tmpreads2.name)

    def assemble(self,):
	
	import glob
	import os
	import sys
	import subprocess
	import multiprocessing

        # find script path
        path = os.path.abspath(sys.argv[0])
        path = os.path.abspath('/'.join(sys.argv[0].split('/')[:-2]))

	try: os.mkdir(self.dataPath+'/assembly')
	except OSError: pass

	sys.stderr.write('assebly of sample '+self.name+'...\n')
	
	
	command = [path+'/bin/velvet_1.2.10/velveth',self.dataPath+'/assembly','21','-fastq','-short',self.dataPath+'/mapqMoreThan20.fq']
	sys.stderr.write('velveth sample '+self.name+'...')# ('+' '.join(command)+') ')
	velveth = subprocess.Popen(command,stdout=subprocess.PIPE,stderr=subprocess.PIPE )
	errdata = velveth.communicate()
	if velveth.returncode != 0:
	    print 'cmd: '+' '.join( command )
	    print 'velveth view Error code', velveth.returncode, errdata
	    sys.exit()
	sys.stderr.write('done\n')
	
	command = [path+'/bin/velvet_1.2.10/velvetg',self.dataPath+'/assembly/','-cov_cutoff','auto']
	sys.stderr.write('velvetg sample '+self.name+'...')#  ('+' '.join(command)+') ')
	velvetg = subprocess.Popen(command,stdout=subprocess.PIPE,stderr=subprocess.PIPE )
	errdata = velvetg.communicate()
	if velvetg.returncode != 0:
	    print 'cmd: '+' '.join( command )
	    print 'velvetg view Error code', velvetg.returncode, errdata
	    sys.exit()
	sys.stderr.write('done\n')
	
	sys.stderr.write('assebly of sample '+self.name+' finished.\n')

    def report(self,):
#	cat $i/contigs.fa
#	mv -v $i* $i
#	echo -e "\n### The reads mapped to the following gi: "
#	samtools view -F 4 -q 20 $i/$i.sorted.bam | cut -f3 | sort | uniq -c;
#	echo -e "\n### same as to the following sequences: "
#	echo -e $(for name in $(samtools view -F 4 -q 20 $i/$i.sorted.bam | cut -f3 | sort | uniq); do grep $name references/ncbi.T_cell_receptor_beta.Human.mRNA.fa;done) | sed s/\>/\\\n/g | grep -P "T|cell|receptor|beta|$" --color=Always
	pass
    

if __name__ == "__main__": main()