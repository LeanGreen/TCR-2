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

    def trimFastqs(self):
	for filePairId,readCount,fastq1,fastq2,sampleId in self.getFastqs():

	    import sys
	    import time
	    try: project = sys.argv[3]
	    except IndexError:
		msg = '#ERROR_MSG#'+time.strftime("%Y-%m-%d:%H:%M:%S",time.localtime())+'# You must give a project number for the creation of sbatch scripts, exiting.\n'
		AnalysisPipe.logfile.write(msg)
		sys.stderr.write(msg)
		sys.exit(1)

	    #
	    # sbatch header
	    #
	    output = ''
	    output += '#! /bin/bash -l'+'\n'
	    output += '#SBATCH -A '+project+'\n'
	    output += '#SBATCH -n 2 -p core'+'\n'
	    output += '#SBATCH -t 72:00:00'+'\n'
	    output += '#SBATCH -J trim.'+self.name+'.'+str(filePairId)+'\n'
	    output += '#SBATCH -e '+self.logPath+'/stderr.trimming.'+self.name+'.'+str(filePairId)+'.txt'+'\n'
	    output += '#SBATCH -o '+self.logPath+'/stdout.trimming.'+self.name+'.'+str(filePairId)+'.txt'+'\n'
	    
	    try:
		output += '#SBATCH --mail-type=All'+'\n'
		output += '#SBATCH --mail-user='+sys.argv[5]+'\n'
	    except IndexError: pass
	    
	    #
	    # define variebles and go to path
	    #
	    output += 'echo "$(date) Running on: $(hostname)"'+'\n'
	    output += 'cd '+self.path+'\n'
	    output += 'echo "-----"'+'\n'
	    
	    output += 'module load bioinfo-tools FastQC #cutadapt/1.5.0'+'\n'
	    output += 'workon py2.7\n\n'

	    #
	    # WGA adapter trimming
	    #
	    output += ''+AnalysisPipe.scriptPath+'/wgaAdapterTrimmer.py -i '+fastq1+' > '+self.dataPath+'/'+str(filePairId)+'.r1.wgaTrimmed.fq 2> '+self.logPath+'/rubiconWgaTrimming.'+str(filePairId)+'.r1.log.txt &\n'
	    output += ''+AnalysisPipe.scriptPath+'/wgaAdapterTrimmer.py -i '+fastq2+' > '+self.dataPath+'/'+str(filePairId)+'.r2.wgaTrimmed.fq 2> '+self.logPath+'/rubiconWgaTrimming.'+str(filePairId)+'.r2.log.txt &\n'
	    output += 'wait\n'
	    
	    output += '\n'
	    output += 'cutadapt -n 10 -g GTGAGTGATGGTTGAGGTAGTGTGGAG -a CTCCACACTACCTCAACCATCACTCAC '+self.dataPath+'/'+str(filePairId)+'.r1.wgaTrimmed.fq > '+self.dataPath+'/'+str(filePairId)+'.r1.wgaTrimmed2.fq  2> '+self.logPath+'/malbacWgaTrimming.'+str(filePairId)+'.r1.log.txt &\n'
	    output += 'cutadapt -n 10 -g GTGAGTGATGGTTGAGGTAGTGTGGAG -a CTCCACACTACCTCAACCATCACTCAC '+self.dataPath+'/'+str(filePairId)+'.r2.wgaTrimmed.fq > '+self.dataPath+'/'+str(filePairId)+'.r2.wgaTrimmed2.fq  2> '+self.logPath+'/malbacWgaTrimming.'+str(filePairId)+'.r2.log.txt &\n'
	    output += 'wait\n'
	    output += '\n'
	    #output += 'rm -v '+self.dataPath+'/'+str(filePairId)+'.r1.wgaTrimmed.fq '+self.dataPath+'/'+str(filePairId)+'.r2.wgaTrimmed.fq\n'
	    
	    #
	    # illumina  adapter trimming
	    #
	    adaptersToTrim = '-a CTGTCTCTTATACACATCTGACGCTGCCGACGA -a CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT'
	    output += 'cutadapt -n 3 '+adaptersToTrim+' '+self.dataPath+'/'+str(filePairId)+'.r1.wgaTrimmed2.fq > '+self.dataPath+'/'+str(filePairId)+'.r1.wgaAndilluminaTrimmed.fq 2> '+self.logPath+'/illuminaAndNexteraTrimming.'+str(filePairId)+'.r1.log.txt &\n'
	    output += 'cutadapt -n 3 '+adaptersToTrim+' '+self.dataPath+'/'+str(filePairId)+'.r2.wgaTrimmed2.fq > '+self.dataPath+'/'+str(filePairId)+'.r2.wgaAndilluminaTrimmed.fq 2> '+self.logPath+'/illuminaAndNexteraTrimming.'+str(filePairId)+'.r2.log.txt &\n'
	    output += 'wait\n'

	    #
	    # remove temp files
	    #
	    #output += 'rm -v '+self.dataPath+'/'+str(filePairId)+'.r1.wgaTrimmed2.fq '+self.dataPath+'/'+str(filePairId)+'.r2.wgaTrimmed2.fq \n'
	    output += 'wait\n'
	    
	    #
	    # quality trimmming
	    #
	    output += ''+AnalysisPipe.scriptPath+'/TrimBWAstyle.pl -q 20 '+self.dataPath+'/'+str(filePairId)+'.r1.wgaAndilluminaTrimmed.fq > '+self.dataPath+'/'+str(filePairId)+'.r1.wgaIlluminaAndQualityTrimmed.fq 2> '+self.logPath+'/qualityTrimming.'+str(filePairId)+'.r1.log.txt &\n'
	    output += ''+AnalysisPipe.scriptPath+'/TrimBWAstyle.pl -q 20 '+self.dataPath+'/'+str(filePairId)+'.r2.wgaAndilluminaTrimmed.fq > '+self.dataPath+'/'+str(filePairId)+'.r2.wgaIlluminaAndQualityTrimmed.fq 2> '+self.logPath+'/qualityTrimming.'+str(filePairId)+'.r2.log.txt &\n'
	    output += 'wait\n'
	    
	    #
	    # remove temp files
	    #
	    #output += 'rm -v '+self.dataPath+'/'+str(filePairId)+'.r1.wgaAndilluminaTrimmed.fq '+self.dataPath+'/'+str(filePairId)+'.r2.wgaAndilluminaTrimmed.fq\n'
	    output += 'wait\n'
	    
	    #
	    # remove empty or "N" only sequences
	    #
	    output += 'python '+AnalysisPipe.scriptPath+'/removeEmptyReads.py '
	    output += self.dataPath+'/'+str(filePairId)+'.r1.wgaIlluminaAndQualityTrimmed.fq '
	    output += self.dataPath+'/'+str(filePairId)+'.r2.wgaIlluminaAndQualityTrimmed.fq '
	    output += self.dataPath+'/'+str(filePairId)+'.r1.allTrimmed.fq '
	    output += self.dataPath+'/'+str(filePairId)+'.r2.allTrimmed.fq '
	    output += self.dataPath+'/'+str(filePairId)+'.singletts.fq '
	    output += '>&2 2> '+self.logPath+'/removeEmptyReads.'+str(filePairId)+'.log.txt\n'
	    
	    #
	    # remove temp files
	    #
	    #output += 'rm -v '+self.dataPath+'/'+str(filePairId)+'.r1.wgaIlluminaAndQualityTrimmed.fq '+self.dataPath+'/'+str(filePairId)+'.r2.wgaIlluminaAndQualityTrimmed.fq\n'
	    output += 'wait\n'
	    
	    #
	    # compress files
	    #
	    output += 'gzip -v9 '+self.dataPath+'/'+str(filePairId)+'.r1.allTrimmed.fq &\n'
	    output += 'gzip -v9 '+self.dataPath+'/'+str(filePairId)+'.r2.allTrimmed.fq  &\n'
	    output += 'gzip -v9 '+self.dataPath+'/'+str(filePairId)+'.singletts.fq &\n'
	    output += 'wait\n'
	    
	    #
	    # FASTQC
	    #
	    output += 'fastqc '+self.dataPath+'/'+str(filePairId)+'.r1.allTrimmed.fq.gz &\n'
	    output += 'fastqc '+self.dataPath+'/'+str(filePairId)+'.r2.allTrimmed.fq.gz &\n'
	    output += 'fastqc '+self.dataPath+'/'+str(filePairId)+'.singletts.fq.gz &\n'
	    output += 'wait\n'
	    output += 'mv -v '+self.dataPath+'/*fastqc* '+self.fastqcPath+'/\n'
	    
	    #
	    # Final output and write script to file
	    #
            output += '\n'+AnalysisPipe.programPath+' '+AnalysisPipe.path+' report\n'
	    output += 'echo'+'\n'
	    output += 'wait'+'\n'
	    output += 'echo "$(date) AllDone"'+'\n'
	    output += 'echo "$(date) AllDone" >&2'+'\n'
	    with open(self.scriptPath+'/trimming.'+self.name+'.'+str(filePairId)+'.sh','w') as outfile: outfile.write(output)

    def mapFastqs(self):

	for filePairId,readCount,fastq1,fastq2,sampleId in self.getFastqs():

	    import sys
	    import time
	    try: project = sys.argv[3]
	    except IndexError:
		msg = '#ERROR_MSG#'+time.strftime("%Y-%m-%d:%H:%M:%S",time.localtime())+'# You must give a project number for the creation of sbatch scripts, exiting.\n'
		AnalysisPipe.logfile.write(msg)
		sys.stderr.write(msg)
		sys.exit(1)

	    #
	    # sbatch header
	    #
	    output = '#! /bin/bash -l'+'\n'
	    output += '#SBATCH -A '+project+'\n'
	    output += '#SBATCH -n 16 -p node'+'\n'
	    output += '#SBATCH -t 72:00:00'+'\n'
	    output += '#SBATCH -J map.'+self.name+'.'+str(filePairId)+'\n'
	    output += '#SBATCH -e '+self.logPath+'/stderr.mapping.'+self.name+'.'+str(filePairId)+'.txt'+'\n'
	    output += '#SBATCH -o '+self.logPath+'/stdout.mapping.'+self.name+'.'+str(filePairId)+'.txt'+'\n'

	    try:
		output += '#SBATCH --mail-type=All'+'\n'
		output += '#SBATCH --mail-user='+sys.argv[5]+'\n'
	    except IndexError: pass

	    #
	    # define variebles and go to path
	    #
	    output += 'echo "$(date) Running on: $(hostname)"'+'\n'
	    output += 'cd '+self.path+'\n'
	    output += 'echo'+'\n'

	    #
	    # Bowtie2 mapping
	    #output += 'module load bioinfo-tools bwa/0.7.8\n'
	    #output += 'bwa mem -t 16 /sw/data/uppnex/reference/Homo_sapiens/GRCh37/program_files/bwa/concat.fa '+self.r1files[0]+' '+self.r2files[0]+' > '+self.sam+'\n'
	    #output += 'bowtie2 -1 '+self.r1files[0]+' -2 '+self.r2files[0]+' --very-sensitive-local -p16 -x '+self.reference+' > '+self.sam+'\n'
	    output += 'bowtie2 --maxins 2000 -p16 '
	    output += '-1 '+self.dataPath+'/'+str(filePairId)+'.r1.allTrimmed.fq.gz '
	    output += '-2 '+self.dataPath+'/'+str(filePairId)+'.r2.allTrimmed.fq.gz '
	    output += '-x '+AnalysisPipe.bowtie2Reference+' '
	    output += '> '+self.dataPath+'/'+str(filePairId)+'.sam '
	    output += '2> '+self.logPath+'/stderr.bowtie2.'+str(filePairId)+'.txt \n'
	    output += 'echo -e "mapping Done. $(date) Running on: $(hostname)" 1>&2'+'\n'
	    #output += 'rm -v '+self.dataPath+'/'+str(filePairId)+'.r1.allTrimmed.fq.gz'+'\n'
	    #output += 'rm -v '+self.dataPath+'/'+str(filePairId)+'.r2.allTrimmed.fq.gz'+'\n'
	    

	    #
	    # Final output and write script to file
	    #
	    output += '\n'+AnalysisPipe.programPath+' '+AnalysisPipe.path+' report\n'
	    output += 'wait'+'\n'
	    output += 'echo "$(date) AllDone"'+'\n'
	    output += 'echo "$(date) AllDone" >&2'+'\n'
	    with open(self.scriptPath+'/mapping.'+self.name+'.'+str(filePairId)+'.sh','w') as outfile: outfile.write(output)

    def mergeMapped(self):
	import sys
	import time
	try: project = sys.argv[3]
	except IndexError:
	    msg = '#ERROR_MSG#'+time.strftime("%Y-%m-%d:%H:%M:%S",time.localtime())+'# You must give a project number for the creation of sbatch scripts, exiting.\n'
	    AnalysisPipe.logfile.write(msg)
	    sys.stderr.write(msg)
	    sys.exit(1)
	try: wgsOrExome = sys.argv[4]
	except IndexError:
	    msg = '#ERROR_MSG#'+time.strftime("%Y-%m-%d:%H:%M:%S",time.localtime())+'# You must give a type for the anaylis wgs or exome, now exiting.\n'
	    AnalysisPipe.logfile.write(msg)
	    sys.stderr.write(msg)
	    sys.exit(1)
	#
	# sbatch header
	#
	output = '#! /bin/bash -l'+'\n'
	output += '#SBATCH -A '+project+'\n'
	output += '#SBATCH -n 1 -p core'+'\n'
	output += '#SBATCH -t 5:00:00'+'\n'
	output += '#SBATCH -J merge.'+self.name+'\n'
	output += '#SBATCH -e '+self.logPath+'/stderr.merge.'+self.name+'.txt'+'\n'
	output += '#SBATCH -o '+self.logPath+'/stdout.merge.'+self.name+'.txt'+'\n'

	try:
	    output += '#SBATCH --mail-type=All'+'\n'
	    output += '#SBATCH --mail-user='+sys.argv[5]+'\n'
	except IndexError: pass

	#
	# define variebles and go to path
	#
	output += 'echo "$(date) Running on: $(hostname)"'+'\n'
	output += 'cd '+self.path+'\n'
	output += 'echo'+'\n'

	#
	# merge
	#
	inputFiles = ' INPUT='+' INPUT='.join([self.dataPath+'/'+str(filePairId)+'.sam' for filePairId,readCount,fastq1,fastq2,sampleId in self.getFastqs()])
	output += 'java -Xmx5g -jar '+AnalysisPipe.picardLocation+'/MergeSamFiles.jar '+inputFiles+' OUTPUT='+self.dataPath+'/'+self.name+'.merged.sam '
	output += '1>&2  2>  '+self.logPath+'/stderr.merging.'+self.name+'.txt \n'
	output += 'echo -e "mapping Done. $(date) Running on: $(hostname)" 1>&2'+'\n'

	#
	# Final output and write script to file
	#
	output += 'wait'+'\n'
	output += 'echo "$(date) AllDone"'+'\n'
	output += 'echo "$(date) AllDone" >&2'+'\n'
	with open(self.scriptPath+'/mergeMapped.'+self.name+'.sh','w') as outfile: outfile.write(output)

    def filterAndFixMerged(self):
	import sys
	import time
	try: project = sys.argv[3]
	except IndexError:
	    msg = '#ERROR_MSG#'+time.strftime("%Y-%m-%d:%H:%M:%S",time.localtime())+'# You must give a project number for the creation of sbatch scripts, exiting.\n'
	    AnalysisPipe.logfile.write(msg)
	    sys.stderr.write(msg)
	    sys.exit(1)
	try: wgsOrExome = sys.argv[4]
	except IndexError:
	    msg = '#ERROR_MSG#'+time.strftime("%Y-%m-%d:%H:%M:%S",time.localtime())+'# You must give a type for the anaylis wgs or exome, now exiting.\n'
	    AnalysisPipe.logfile.write(msg)
	    sys.stderr.write(msg)
	    sys.exit(1)

        #
        # sbatch header
        #
	output = '#! /bin/bash -l'+'\n'
	output += '#SBATCH -A '+project+'\n'
	output += '#SBATCH -n 1 -p core'+'\n'
	output += '#SBATCH -t 72:00:00'+'\n'
	output += '#SBATCH -J fnf.'+self.name+'\n'
	output += '#SBATCH -e '+self.logPath+'/stderr.filterAndFix.'+self.name+'.txt'+'\n'
	output += '#SBATCH -o '+self.logPath+'/stdout.filterAndFix.'+self.name+'.txt'+'\n'

	try:
	    output += '#SBATCH --mail-type=All'+'\n'
	    output += '#SBATCH --mail-user='+sys.argv[5]+'\n'
	except IndexError: pass

        #
        # define variebles and go to path
        #
        output += 'echo "$(date) Running on: $(hostname)"'+'\n'
        output += 'cd '+self.path+'\n'
        output += 'echo'+'\n'
        
        #
        # convert to bam file
        #
        output += 'java -Xmx5g -jar '+AnalysisPipe.picardLocation+'/SamFormatConverter.jar MAX_RECORDS_IN_RAM=2500000 '
	output += 'INPUT='+ self.dataPath+'/'+self.name+'.merged.sam '
	output += 'OUTPUT='+self.dataPath+'/'+self.name+'.merged.bam '
	output += '1>&2  2> '+self.logPath+'/stderr.sam2bam.'+self.name+'.txt \n'
        output += 'echo -e "sam2bam Done. $(date) Running on: $(hostname)" 1>&2'+'\n'
        #output += 'rm -v '+self.dataPath+'/'+self.name+'.merged.sam\n'

        #
        # sort the bam file
        #
        output += 'java -Xmx5g -jar '+AnalysisPipe.picardLocation+'/SortSam.jar MAX_RECORDS_IN_RAM=2500000 SORT_ORDER=coordinate '
	output += 'INPUT='+ self.dataPath+'/'+self.name+'.merged.bam '
        output += 'OUTPUT='+self.dataPath+'/'+self.name+'.sorted.bam '
	output += 'CREATE_INDEX=true 1>&2  2> '
	output += self.logPath+'/stderr.sortBam.'+self.name+'.txt \n'
        output += 'echo -e "bam2sort Done. $(date) Running on: $(hostname)" 1>&2'+'\n'
        #output += 'rm -v '+self.dataPath+'/'+self.name+'.merged.bam\n'

        #
        # mark duplicates
        #
        output += 'java -Xmx5g -jar '+AnalysisPipe.picardLocation+'/MarkDuplicates.jar MAX_RECORDS_IN_RAM=2500000 VALIDATION_STRINGENCY=LENIENT '
	output += 'INPUT='+ self.dataPath+'/'+self.name+'.sorted.bam '
	output += 'OUTPUT='+self.dataPath+'/'+self.name+'.marked.bam '
	output += 'METRICS_FILE='+self.logPath+'/markDuplicatesMetrix.'+self.name+'.txt '
	output += '1>&2  2> '+self.logPath+'/stderr.markDuplicates.'+self.name+'.txt \n'
        output += 'echo -e "mark Done. $(date) Running on: $(hostname)" 1>&2'+'\n'
        #output += 'rm -v '+self.dataPath+'/'+self.name+'.sorted.bam\n'

        #
        # fix missing information
        #
        output += 'java -Xmx5g -jar '+AnalysisPipe.picardLocation+'/AddOrReplaceReadGroups.jar '
        output += 'MAX_RECORDS_IN_RAM=2500000 '
        output += 'INPUT='+ self.dataPath+'/'+self.name+'.marked.bam '
        output += 'OUTPUT='+self.dataPath+'/'+self.name+'.fixed.bam '
        output += 'CREATE_INDEX=true RGID='+self.name+' RGLB='+self.name+' RGPL=ILLUMINA RGSM='+self.name+' RGCN="NA" RGPU="NA"'+'  '
	output += '1>&2  2> '+self.logPath+'/stderr.addAndReplaceReadGroups.'+self.name+'.txt \n'
        output += 'echo "addorreplace Done. $(date) Running on: $(hostname)" 1>&2'+'\n'
        #output += 'rm -v '+self.dataPath+'/'+self.name+'.marked.bam\n'

        #
        # samtools flagstat
        #
        output += 'samtools flagstat '+self.dataPath+'/'+self.name+'.fixed.bam'+' > '+self.logPath+'/fixedBamFlagstat.'+self.name+'.txt \n'
        output += 'echo "flagstat Done. $(date) Running on: $(hostname)" 1>&2'+'\n'

        #
        # Final output and write script to file
        #
        output += '\n'+AnalysisPipe.programPath+' '+AnalysisPipe.path+' report\n'
        output += 'wait'+'\n'
        output += 'echo "$(date) AllDone"'+'\n'
        output += 'echo "$(date) AllDone" >&2'+'\n'
	with open(self.scriptPath+'/FilterAndFix.'+self.name+'.sh','w') as outfile: outfile.write(output)

    def realignerTargetCreator(self):
	import sys
	import time
	try: project = sys.argv[3]
	except IndexError:
	    msg = '#ERROR_MSG#'+time.strftime("%Y-%m-%d:%H:%M:%S",time.localtime())+'# You must give a project number for the creation of sbatch scripts, exiting.\n'
	    AnalysisPipe.logfile.write(msg)
	    sys.stderr.write(msg)
	    sys.exit(1)
	try: wgsOrExome = sys.argv[4]
	except IndexError:
	    msg = '#ERROR_MSG#'+time.strftime("%Y-%m-%d:%H:%M:%S",time.localtime())+'# You must give a type for the anaylis wgs or exome, now exiting.\n'
	    AnalysisPipe.logfile.write(msg)
	    sys.stderr.write(msg)
	    sys.exit(1)

        #
        # sbatch header
        #

	output = '#! /bin/bash -l'+'\n'
	output += '#SBATCH -A '+project+'\n'
	output += '#SBATCH -n 16 -p node'+'\n'
	output += '#SBATCH -t 72:00:00'+'\n'
	output += '#SBATCH -J realTC.'+self.name+'\n'
	output += '#SBATCH -e '+self.logPath+'/stderr.realTC.'+self.name+'.txt'+'\n'
	output += '#SBATCH -o '+self.logPath+'/stdout.realTC.'+self.name+'.txt'+'\n'

	try:
	    output += '#SBATCH --mail-type=All'+'\n'
	    output += '#SBATCH --mail-user='+sys.argv[5]+'\n'
	except IndexError: pass

        #
        # define variebles and go to path
        #
        output += 'echo "$(date) Running on: $(hostname)"\n'

        #
        # Find targets for indel realignment
        #
        output += 'echo -e "-> RealignerTargetCreator <-"\n'
        output += 'java -Xmx72g -jar '+AnalysisPipe.gatkLocation+' -T RealignerTargetCreator '
	output += '-nt 16 '
	output += '-I '+self.dataPath+'/'+self.name+'.fixed.bam'+' '
	output += '-R '+AnalysisPipe.bowtie2Reference+' '
	output += '-o '+self.dataPath+'/'+self.name+'.reAlignemntTargetIntervals.bed '
        output += ' -known '+AnalysisPipe.gatkBundleLocation+'/Mills_and_1000G_gold_standard.indels.b37.vcf'
        output += ' -known '+AnalysisPipe.gatkBundleLocation+'/1000G_phase1.indels.b37.vcf '
	output += '1>&2 2> '+self.logPath+'/stderr.RealignerTargetCreator.'+self.name+'.txt;'
        output += '\n'
        
        with open(self.scriptPath+'/realignerTargetCreator.'+self.name+'.sh','w') as outfile: outfile.write(output)
    
    def reAlignAndReCalibrate(self):
	import sys
	import time
	try: project = sys.argv[3]
	except IndexError:
	    msg = '#ERROR_MSG#'+time.strftime("%Y-%m-%d:%H:%M:%S",time.localtime())+'# You must give a project number for the creation of sbatch scripts, exiting.\n'
	    AnalysisPipe.logfile.write(msg)
	    sys.stderr.write(msg)
	    sys.exit(1)
	try: wgsOrExome = sys.argv[4]
	except IndexError:
	    msg = '#ERROR_MSG#'+time.strftime("%Y-%m-%d:%H:%M:%S",time.localtime())+'# You must give a type for the anaylis wgs or exome, now exiting.\n'
	    AnalysisPipe.logfile.write(msg)
	    sys.stderr.write(msg)
	    sys.exit(1)

        #
        # sbatch header
        #
	output = '#! /bin/bash -l'+'\n'
	output += '#SBATCH -A '+project+'\n'
	output += '#SBATCH -n 1 -p core'+'\n'
	output += '#SBATCH -t 72:00:00'+'\n'
	output += '#SBATCH -J reAlign.'+self.name+'\n'
	output += '#SBATCH -e '+self.logPath+'/stderr.reAlign.'+self.name+'.txt'+'\n'
	output += '#SBATCH -o '+self.logPath+'/stdout.reAlign.'+self.name+'.txt'+'\n'

	try:
	    output += '#SBATCH --mail-type=All'+'\n'
	    output += '#SBATCH --mail-user='+sys.argv[5]+'\n'
	except IndexError: pass
        output += 'echo "$(date) Running on: $(hostname)"\n'

        #
        # Realign reads around indels
        #
        output += 'echo -e "-> IndelRealigner <-"\n'
        output += 'java -Xmx5g -jar '+AnalysisPipe.gatkLocation+' -T IndelRealigner '
	output += '-I '+self.dataPath+'/'+self.name+'.fixed.bam'+' '
	output += '-R '+AnalysisPipe.bowtie2Reference+' '
	output += '-targetIntervals '+self.dataPath+'/'+self.name+'.reAlignemntTargetIntervals.bed '
        output += ' -o '+self.dataPath+'/'+self.name+'.reAligned.bam'+' '
        output += ' -known '+AnalysisPipe.gatkBundleLocation+'/Mills_and_1000G_gold_standard.indels.b37.vcf'
        output += ' -known '+AnalysisPipe.gatkBundleLocation+'/1000G_phase1.indels.b37.vcf  '
	output += '1>&2 2> '+self.logPath+'/stderr.indelRealigner.'+self.name+'.txt;'+'\n'
        output += '\n'
        output += 'echo "Done. $(date) Running on: $(hostname)"\n'
        output += 'echo "$(date) Running on: $(hostname)"\n'
        #output += 'rm -v '+self.dataPath+'/'+self.name+'.fixed.bam'+'\n'
        
        #
        # Quality recalibration
        #
        output += 'echo -e "-> BaseRecalibrator <-"\n'
        output += 'java -Xmx5g -jar '+AnalysisPipe.gatkLocation+' -T BaseRecalibrator '
	output += '-I '+self.dataPath+'/'+self.name+'.reAligned.bam'+' '
	output += '-R '+AnalysisPipe.bowtie2Reference+' '
	output += '-o '+self.dataPath+'/'+self.name+'.BQSR.grp'+' '
        output += ' -knownSites '+AnalysisPipe.gatkBundleLocation+'/dbsnp_138.b37.vcf '
	output += '1>&2 2> '+self.logPath+'/stderr.baseRecalibrator.'+self.name+'.txt;'+'\n'

        output += '\n'
        output += 'echo -e "-> PrintReads <-"\n'
        output += 'java -Xmx5g -jar '+AnalysisPipe.gatkLocation+' -T PrintReads '
	output += '-I '+self.dataPath+'/'+self.name+'.reAligned.bam'+' '
	output += '-R '+AnalysisPipe.bowtie2Reference+' '
	output += '-BQSR '+self.dataPath+'/'+self.name+'.BQSR.grp'+' '
	output += '-o '+self.dataPath+'/'+self.name+'.reCalibrated.bam'+' '
	output += '1>&2 2> '+self.logPath+'/stderr.printreads.txt ;\n'
        #output += 'rm -v '+self.dataPath+'/'+self.name+'.reAligned.bam'+'\n'
        output += 'samtools flagstat '+self.dataPath+'/'+self.name+'.reCalibrated.bam > '+self.logPath+'/reCalibratedBamFlagstat.'+self.name+'.txt \n'

	output += 'samtools view -b -F 4 '   +self.dataPath+'/'+self.name+'.reCalibrated.bam > '+self.dataPath+'/'+self.name+'.unmapRemoved.bam  2> '+self.logPath+'/stderr.samtoolsView.removeUnmap.'+self.name+'.txt \n'
	output += 'java -Xmx5g -jar '+AnalysisPipe.picardLocation+'/BuildBamIndex.jar INPUT='+self.dataPath+'/'+self.name+'.unmapRemoved.bam '+'1>&2  2>  '+self.logPath+'/stderr.buildIndex1.'+self.name+'.txt \n'
	output += 'samtools flagstat '+self.dataPath+'/'+self.name+'.unmapRemoved.bam > '+self.logPath+'/unmapRemovedBamFlagstat.'+self.name+'.txt \n'
        #output += 'rm -v '+self.dataPath+'/'+self.name+'.reCalibrated.bam\n'

	output += 'samtools view -b -q 20 '  +self.dataPath+'/'+self.name+'.unmapRemoved.bam > '+self.dataPath+'/'+self.name+'.qualFiltered.bam  2> '+self.logPath+'/stderr.samtoolsView.qualFilter.'+self.name+'.txt \n'
	output += 'java -Xmx5g -jar '+AnalysisPipe.picardLocation+'/BuildBamIndex.jar INPUT='+self.dataPath+'/'+self.name+'.qualFiltered.bam '+'1>&2  2>  '+self.logPath+'/stderr.buildIndex2.'+self.name+'.txt \n'
	output += 'samtools flagstat '+self.dataPath+'/'+self.name+'.qualFiltered.bam > '+self.logPath+'/qualFilteredBamFlagstat.'+self.name+'.txt \n'
	#output += 'rm -v '+self.dataPath+'/'+self.name+'.unmapRemoved.bam\n'

	output += 'samtools view -b -F 1024 '+self.dataPath+'/'+self.name+'.qualFiltered.bam > '+self.dataPath+'/'+self.name+'.noDuplicates.bam  2> '+self.logPath+'/stderr.samtoolsView.removeDups.'+self.name+'.txt \n'
	output += 'java -Xmx5g -jar '+AnalysisPipe.picardLocation+'/BuildBamIndex.jar INPUT='+self.dataPath+'/'+self.name+'.noDuplicates.bam '+'1>&2  2>  '+self.logPath+'/stderr.buildIndex3.'+self.name+'.txt \n'
	output += 'samtools flagstat '+self.dataPath+'/'+self.name+'.noDuplicates.bam > '+self.logPath+'/noDuplicatesBamFlagstat.'+self.name+'.txt \n'
	#output += 'rm -v '+self.dataPath+'/'+self.name+'.qualFiltered.bam\n'

        #
        # Final output and write script to file
        #
        output += '\n'+AnalysisPipe.programPath+' '+AnalysisPipe.path+' report\n'
        output += 'echo "Done. $(date) Running on: $(hostname)"\n'
        output += 'wait\n'
        output += 'echo "$(date) AllDone"\n'
        with open(self.scriptPath+'/reAlignAndReCalibrate.'+self.name+'.sh','w') as outfile: outfile.write(output)

    def haplotypeCalling(self):
	import sys
	import time
	try: project = sys.argv[3]
	except IndexError:
	    msg = '#ERROR_MSG#'+time.strftime("%Y-%m-%d:%H:%M:%S",time.localtime())+'# You must give a project number for the creation of sbatch scripts, exiting.\n'
	    AnalysisPipe.logfile.write(msg)
	    sys.stderr.write(msg)
	    sys.exit(1)
	try: wgsOrExome = sys.argv[4]
	except IndexError:
	    msg = '#ERROR_MSG#'+time.strftime("%Y-%m-%d:%H:%M:%S",time.localtime())+'# You must give a type for the anaylis wgs or exome, now exiting.\n'
	    AnalysisPipe.logfile.write(msg)
	    sys.stderr.write(msg)
	    sys.exit(1)

        #
        # sbatch header
        #
	output = '#! /bin/bash -l'+'\n'
	output += '#SBATCH -A '+project+'\n'
	output += '#SBATCH -n 1 -p core'+'\n'
	output += '#SBATCH -t 72:00:00'+'\n'
	output += '#SBATCH -J hapCal.'+self.name+'\n'
	output += '#SBATCH -e '+self.logPath+'/stderr.haplotypeCalling.'+self.name+'.txt'+'\n'
	output += '#SBATCH -o '+self.logPath+'/stdout.haplotypeCalling.'+self.name+'.txt'+'\n'

	try:
	    output += '#SBATCH --mail-type=All'+'\n'
	    output += '#SBATCH --mail-user='+sys.argv[5]+'\n'
	except IndexError: pass
        output += 'echo "$(date) Running on: $(hostname)"\n'

        output += 'echo "HC" '+'\n'
        
        output += 'java -Xmx5g -jar '+AnalysisPipe.gatkLocation+' '
        output += '-T HaplotypeCaller '
	output += '-R '+AnalysisPipe.bowtie2Reference+' '
        output += '-I '+self.dataPath+'/'+self.name+'.noDuplicates.bam '
        output += '--genotyping_mode DISCOVERY '
        output += '-stand_emit_conf 10 '
        output += '-stand_call_conf 30 '
        if wgsOrExome == 'exome': output += '-L '+AnalysisPipe.referencePath+'/truseq_exome_targeted_regions.hg19.bed.chr.columnReOrdered.withHeader.chrRem.bed '
	elif wgsOrExome == 'wgs': output += '-L '+AnalysisPipe.referencePath+'/wgs.bed '
        output += '--dbsnp '+AnalysisPipe.gatkBundleLocation+'/dbsnp_138.b37.vcf '
        output += '--annotation AlleleBalance --annotation AlleleBalanceBySample --annotation BaseCounts --annotation BaseQualityRankSumTest '
        output += '--annotation ChromosomeCounts --annotation ClippingRankSumTest --annotation Coverage --annotation DepthPerAlleleBySample '
        output += '--annotation DepthPerSampleHC --annotation FisherStrand --annotation GCContent --annotation HaplotypeScore --annotation HardyWeinberg '
        output += '--annotation HomopolymerRun --annotation InbreedingCoeff --annotation LikelihoodRankSumTest --annotation LowMQ '
        output += '--annotation MVLikelihoodRatio --annotation MappingQualityRankSumTest --annotation MappingQualityZero --annotation MappingQualityZeroBySample '
        output += '--annotation NBaseCount --annotation QualByDepth --annotation RMSMappingQuality --annotation ReadPosRankSumTest --annotation SampleList '
        output += '--annotation SnpEff --annotation SpanningDeletions --annotation StrandBiasBySample --annotation TandemRepeatAnnotator '
        output += '--annotation TransmissionDisequilibriumTest --annotation VariantType '#--annotation StrandOddsRatio 
        output += '--emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 '
        output += '-o '+self.dataPath+'/'+self.name+'.gvcf '
	output += '1>&2 2> '+self.logPath+'/stderr.haplotypeCallerGatk.'+self.name+'.txt &'+'\n'
        output += 'wait'+'\n'

        #
        # Final output and write script to file
        #
        output += 'echo "Done. $(date) Running on: $(hostname)"\n'
        output += 'wait\n'
        output += 'echo "$(date) AllDone"\n'
        with open(self.scriptPath+'/haplotypeCalling.'+self.name+'.sh','w') as outfile: outfile.write(output)

    def qcSteps(self):
	import sys
	import time
	try: project = sys.argv[3]
	except IndexError:
	    msg = '#ERROR_MSG#'+time.strftime("%Y-%m-%d:%H:%M:%S",time.localtime())+'# You must give a project number for the creation of sbatch scripts, exiting.\n'
	    AnalysisPipe.logfile.write(msg)
	    sys.stderr.write(msg)
	    sys.exit(1)
	try: wgsOrExome = sys.argv[4]
	except IndexError:
	    msg = '#ERROR_MSG#'+time.strftime("%Y-%m-%d:%H:%M:%S",time.localtime())+'# You must give a type for the anaylis wgs or exome, now exiting.\n'
	    AnalysisPipe.logfile.write(msg)
	    sys.stderr.write(msg)
	    sys.exit(1)

        #
        # sbatch header
        #
	output = '#! /bin/bash -l'+'\n'
	output += '#SBATCH -A '+project+'\n'
	output += '#SBATCH -n 1 -p core'+'\n'
	output += '#SBATCH -t 72:00:00'+'\n'
	output += '#SBATCH -J qcSteps.'+self.name+'\n'
	output += '#SBATCH -e '+self.logPath+'/stderr.qcSteps.'+self.name+'.txt'+'\n'
	output += '#SBATCH -o '+self.logPath+'/stdout.qcSteps.'+self.name+'.txt'+'\n'

	try:
	    output += '#SBATCH --mail-type=All'+'\n'
	    output += '#SBATCH --mail-user='+sys.argv[5]+'\n'
	except IndexError: pass
        output += 'echo "$(date) Running on: $(hostname)"\n'

        #
        # GATK callable Loci
        #
        output += 'echo "$(date) Running on: $(hostname)"'+'\n'
        output += 'echo -e "-> CallableLoci <-"'+'\n'
        output += 'java -Xmx5g -jar '+AnalysisPipe.gatkLocation+' -T CallableLoci '
	output +='-I '+self.dataPath+'/'+self.name+'.noDuplicates.bam '
	output +='-summary '+self.dataPath+'/'+self.name+'.callableLociSummary.txt '
	output +='-o '+self.dataPath+'/'+self.name+'.callableLoci.bed '
	output +='-R '+AnalysisPipe.bowtie2Reference+' '+'\n'
        output += 'echo "Done. $(date) Running on: $(hostname)"'+'\n'
        output += 'echo'+'\n'
        output += 'echo "-----"'+'\n'

        #
        # qacompute
        #
        output += 'echo "$(date) Running on: $(hostname)"'+'\n'
        output += 'echo -e "-> Pauls qacompute <-"'+'\n'
        output += '/proj/b2010052/scripts/qaCompute -d -q 10 '
	output += '-m '+self.dataPath+'/'+self.name+'.noDuplicates.bam '
	output += self.dataPath+'/'+self.name+'.qacompute.out '
	output += '> '+self.logPath+'/'+self.name+'.qacompute.stdout.txt '
	output += '2> '+self.logPath+'/'+self.name+'.qacompute.stderr.txt '+'\n'
        output += 'echo "Done. $(date) Running on: $(hostname)"'+'\n'

        #
        # picard HS metrics
        #
        output += 'java -Xmx3g -jar '+AnalysisPipe.picardLocation+'/CalculateHsMetrics.jar '
	if wgsOrExome == 'exome':output += 'BAIT_INTERVALS='  +AnalysisPipe.referencePath+'/truseq_exome_targeted_regions.hg19.bed.chr.columnReOrdered.withHeader.chrRem '
	elif wgsOrExome == 'wgs':output += 'BAIT_INTERVALS='  +AnalysisPipe.referencePath+'/wgs '
	if wgsOrExome == 'exome':output += 'TARGET_INTERVALS='+AnalysisPipe.referencePath+'/truseq_exome_targeted_regions.hg19.bed.chr.columnReOrdered.withHeader.chrRem '
	elif wgsOrExome == 'wgs':output += 'TARGET_INTERVALS='+AnalysisPipe.referencePath+'/wgs '
	output += 'INPUT='+self.dataPath+'/'+self.name+'.noDuplicates.bam '
	output += 'OUTPUT='+self.dataPath+'/'+self.name+'.hs_metrics.summary.txt '
	output += 'PER_TARGET_COVERAGE='+self.dataPath+'/'+self.name+'.hs_metrics.perTargetCoverage.txt '
	output += 'REFERENCE_SEQUENCE='+AnalysisPipe.bowtie2Reference+'  '
	output += '1>&2 2> '+self.logPath+'/'+self.name+'.stderr.caluclateHsmetrics.txt \n'


	#
	# make files for coverage checks
	#
	if wgsOrExome == 'exome':output += "bedtools coverage -abam "+self.dataPath+'/'+self.name+'.noDuplicates.bam'+" -b "+AnalysisPipe.referencePath+"/truseq_exome_targeted_regions.hg19.bed.chr.columnReOrdered.withHeader.chrRem.bed -d > "+self.dataPath+'/'+self.name+'.bedtools.coverage.bed\n'
	if wgsOrExome == 'exome':output += "bedtools coverage -abam "+self.dataPath+'/'+self.name+'.noDuplicates.bam'+" -b "+AnalysisPipe.referencePath+"/truseq_exome_targeted_regions.hg19.bed.chr.columnReOrdered.withHeader.chrRem.bed -hist > "+self.dataPath+'/'+self.name+'.bedtools.exomecoverage.histogram\n'
	output += "bedtools genomecov -ibam "+self.dataPath+'/'+self.name+'.noDuplicates.bam'+" > "+self.dataPath+'/'+self.name+'.bedtools.genomecov.histogram\n'
	output += "bedtools genomecov -bga -ibam "+self.dataPath+'/'+self.name+'.noDuplicates.bam'+" > "+self.dataPath+'/'+self.name+'.bedtools.genomecov.bedgraph\n'

	if wgsOrExome == 'exome':output += "bedtools coverage -split -abam "+self.dataPath+'/'+self.name+'.noDuplicates.bam'+" -b "+AnalysisPipe.referencePath+"/truseq_exome_targeted_regions.hg19.bed.chr.columnReOrdered.withHeader.chrRem.bed -d > "+self.dataPath+'/'+self.name+'.bedtools.coverage.nonPhysical.bed\n'
	if wgsOrExome == 'exome':output += "bedtools coverage -split -abam "+self.dataPath+'/'+self.name+'.noDuplicates.bam'+" -b "+AnalysisPipe.referencePath+"/truseq_exome_targeted_regions.hg19.bed.chr.columnReOrdered.withHeader.chrRem.bed -hist > "+self.dataPath+'/'+self.name+'.bedtools.exomecoverage.nonPhysical.histogram\n'
	output += "bedtools genomecov -split -ibam "+self.dataPath+'/'+self.name+'.noDuplicates.bam'+" > "+self.dataPath+'/'+self.name+'.bedtools.genomecov.nonPhysical.histogram\n'
	output += "bedtools genomecov -split -bga -ibam "+self.dataPath+'/'+self.name+'.noDuplicates.bam'+" > "+self.dataPath+'/'+self.name+'.bedtools.genomecov.nonPhysical.bedgraph\n'


	output += "awk '{print $7}' "+self.dataPath+'/'+self.name+'.bedtools.coverage.bed'+" | sort -n | uniq -c | awk '{print $1\"\\t\"$2}'> "+self.dataPath+'/'+self.name+".coverageDistribution.tsv\n"
	output += "awk '{print $1 \"\\t\" $2+$6-1 \"\\t\" $7}' "+self.dataPath+'/'+self.name+'.bedtools.coverage.bed'+" > "+self.dataPath+'/'+self.name+".depthPerPosition.tsv\n"

	
        #
        # Final output and write script to file
        #
        output += '\n'+AnalysisPipe.programPath+' '+AnalysisPipe.path+' report\n'
	output += 'echo'+'\n'
        output += 'wait'+'\n'
        output += 'echo "$(date) AllDone"'+'\n'
        output += 'echo "$(date) AllDone" >&2'+'\n'
        with open(self.scriptPath+'/qcSteps.'+self.name+'.sh','w') as outfile: outfile.write(output)

    def updateOrgReadCounts(self,):
	AnalysisPipe.database.updateFastqReadCount(self)

    def getReadOrientationStats(self,):

	import time
	import operator
	import shutil
	import os
	import sys
	import pysam
	import time
	orientations = {'PP':0,'FR':0,'FRsp':0,'FF':0,'RF':0,'RR':0,'difChrom':0}
	
	try: uppmax_temp = os.environ["SNIC_TMP"]
	except:
	    uppmax_temp = None
	    print 'Not on uppmax no temporary directory'
    
	if not os.path.exists(self.dataPath+'/'+self.name+'.noDuplicates.bam'): AnalysisPipe.logfile.write('#WARNING#'+time.strftime("%Y-%m-%d:%H:%M:%S",time.localtime())+'# Skipping oriantations stats for sample '+self.name+' the infile has not been created yet...\n'); return orientations
    
	if uppmax_temp:
	    try:os.mkdir(uppmax_temp+'/fnuttglugg_TMP')
	    except OSError:pass
	    if not os.path.exists(uppmax_temp+'/fnuttglugg_TMP'+'/'+self.name+'.noDuplicates.bam'):
		AnalysisPipe.logfile.write('#LOGMSG#'+time.strftime("%Y-%m-%d:%H:%M:%S",time.localtime())+'# Copying '+self.dataPath+'/'+self.name+'.noDuplicates.bam'+' to temp location for faster reading from disk, '+uppmax_temp+'/fnuttglugg_TMP'+'/'+self.name+'.noDuplicates.bam'+' \n')
		shutil.copy(self.dataPath+'/'+self.name+'.noDuplicates.bam',uppmax_temp+'/fnuttglugg_TMP'+'/'+self.name+'.noDuplicates.bam')
	    else:
		print 'WARNING: rerun of '+uppmax_temp+'/fnuttglugg_TMP'+'/'+self.name+'.noDuplicates.bam'+' skipping copy!!'
		AnalysisPipe.logfile.write('#WARNING#'+time.strftime("%Y-%m-%d:%H:%M:%S",time.localtime())+'# WARNING: rerun of '+uppmax_temp+'/fnuttglugg_TMP'+'/'+self.name+'.noDuplicates.bam'+' skipping copy!!\n')
	    bamfileName  = uppmax_temp+'/fnuttglugg_TMP'+'/'+self.name+'.noDuplicates.bam'
	else:bamfileName = self.dataPath+'/'+self.name+'.noDuplicates.bam'
	if uppmax_temp:
	    try:os.mkdir(uppmax_temp+'/fnuttglugg_TMP')
	    except OSError:pass
	    if not os.path.exists(uppmax_temp+'/fnuttglugg_TMP'+'/'+self.name+'.noDuplicates.bai'):
		AnalysisPipe.logfile.write('#LOGMSG#'+time.strftime("%Y-%m-%d:%H:%M:%S",time.localtime())+'# Copying '+self.dataPath+'/'+self.name+'.noDuplicates.bai'+' to temp location for faster reading from disk, '+uppmax_temp+'/fnuttglugg_TMP'+'/'+self.name+'.noDuplicates.bai'+' \n')
		try: shutil.copy(self.dataPath+'/'+self.name+'.noDuplicates.bai',uppmax_temp+'/fnuttglugg_TMP'+'/'+self.name+'.noDuplicates.bai')
		except IOError as e: pass
	    else:
		print 'WARNING: rerun of '+uppmax_temp+'/fnuttglugg_TMP'+'/'+self.name+'.noDuplicates.bai'+' skipping copy!!'
		AnalysisPipe.logfile.write('#WARNING#'+time.strftime("%Y-%m-%d:%H:%M:%S",time.localtime())+'# WARNING: rerun of '+uppmax_temp+'/fnuttglugg_TMP'+'/'+self.name+'.noDuplicates.bai'+' skipping copy!!\n')
	
#	bamfileName = 	self.dataPath+'/'+self.name+'.noDuplicates.bam'
	try: bamfile = pysam.Samfile(bamfileName, "rb")
	except IOError:
	    AnalysisPipe.logfile.write('#WARNING#'+time.strftime("%Y-%m-%d:%H:%M:%S",time.localtime())+'# Skipping oriantation stats for sample '+self.name+' the infile has not been created yet...\n')
	    return orientations
	except ValueError:
	    AnalysisPipe.logfile.write('#WARNING#'+time.strftime("%Y-%m-%d:%H:%M:%S",time.localtime())+'# Skipping oriantation stats for sample '+self.name+' the infile is not finished for processing...\n')
	    return orientations
    
	AnalysisPipe.logfile.write('#LOGMSG#'+time.strftime("%Y-%m-%d:%H:%M:%S",time.localtime())+'# Loading oriantation stats for sample '+self.name+'...\n')
	try:
	    for read in bamfile.fetch():

		orientation = None
		if read.is_paired and not read.is_unmapped and not read.mate_is_unmapped:
		    if read.tid == read.rnext:
			if read.is_proper_pair:
			    orientation = 'PP'
			elif   read.pos < read.pnext:
			    if read.is_reverse: orientation = 'R'
			    else: orientation = 'F'
			    if read.mate_is_reverse: orientation += 'R'
			    else: orientation += 'F'
			elif read.pnext < read.pos:
			    if read.mate_is_reverse: orientation = 'R'
			    else: orientation = 'F'
			    if read.is_reverse: orientation += 'R'
			    else: orientation += 'F'
			elif   read.pos == read.pnext: orientation = 'FRsp'
		    else: orientation = 'difChrom'
		    orientations[orientation]+=1
	    total = sum(orientations.values())

	    output = self.name+'\n'
	    for thingy,plingy in orientations.iteritems():
		output+= thingy+' '+str(percentage(plingy,total))+'\n'
	    print output

	except ValueError as e:
	    if e == 'fetch called on bamfile without index':
		AnalysisPipe.logfile.write('#WARNING#'+time.strftime("%Y-%m-%d:%H:%M:%S",time.localtime())+'# Skipping insert size plot for sample '+self.name+' the bam index is not present...\n')
		sys.stderr.write('#WARNING#'+time.strftime("%Y-%m-%d:%H:%M:%S",time.localtime())+'# Skipping insert size plot for sample '+self.name+' the bam index is not present...\n')
		return orientations
	AnalysisPipe.logfile.write('#LOGMSG#'+time.strftime("%Y-%m-%d:%H:%M:%S",time.localtime())+'# Loaded oriantation stats for sample '+self.name+' joining to main...\n')
	return orientations

    def getStats(self):
	
	import re
	
	self.updateOrgReadCounts()
	
	stats = {}
	stats['orientations'] = self.getReadOrientationStats()
	stats['rubiconWgaTrimming'] = {}
	stats['malbacWgaTrimming']  = {}
	stats['illuminaAndNexteraTrimming']  = {}
	stats['qualityTrimming']  = {}
	stats['removeEmptyReads']  = {}
	stats['bowtie2'] = {}
	
	for filePairId,readCount,fastq1,fastq2,sampleId in self.getFastqs():
	    stats['rubiconWgaTrimming'][filePairId] = {'r1':None,'r2':None}
	    stats['malbacWgaTrimming'][filePairId]  = {'r1':None,'r2':None}
	    stats['illuminaAndNexteraTrimming'][filePairId]  = {'r1':None,'r2':None}
	    stats['qualityTrimming'][filePairId]    = {'r1':None,'r2':None}
	    
	    for read in ['r1','r2']:
		stats['rubiconWgaTrimming'][filePairId][read]         = extractData(infile=self.logPath+'/rubiconWgaTrimming.'+str(filePairId)+'.'+read+'.log.txt',        pattern="Running wgaAdapterTrimmer.py\nProcessed a total of\t(?P<totalReads>\d+)\treads. \(.+\)\nProcessed a total of\t(?P<totalBases>\d+)\tbases \(.+\).\ntrimmed a total of\t(?P<trimmedBases>\d+)\tbases in the start of reads \(.+\).\nwgaAdapterTrimmer.py done exiting ...\n?")
		stats['malbacWgaTrimming'][filePairId][read]          = extractData(infile=self.logPath+'/malbacWgaTrimming.'+str(filePairId)+'.'+read+'.log.txt',         pattern="cutadapt version .+\nCommand line parameters: -n 10 -g GTGAGTGATGGTTGAGGTAGTGTGGAG -a CTCCACACTACCTCAACCATCACTCAC .+\nMaximum error rate\: .+\%\n\s+No. of adapters\: 2\n\s+Processed reads\:\s+(?P<totalReads>\d+)\n\s+Processed bases\:\s+(?P<totalBases>\d+) bp \(.+ Mbp\)\n\s+Trimmed reads\:\s+(?P<trimmedReads>\d+) \(.+\%\)\n\s+Trimmed bases\:\s+(?P<trimmedBases>\d+) bp \(.+ Mbp\) \(.+\% of total\)\n\s+Too short reads\:\s+.+ \(.+\% of processed reads\)\n\s+Too long reads\:\s+.+ \(.+\% of processed reads\)\n\s+Total time\:\s+.+ s\n\s+Time per read\:\s+.+ ms")
		stats['illuminaAndNexteraTrimming'][filePairId][read] = extractData(infile=self.logPath+'/illuminaAndNexteraTrimming.'+str(filePairId)+'.'+read+'.log.txt',pattern="cutadapt version .+\nCommand line parameters: -n .+\nMaximum error rate\: .+\%\n\s+No. of adapters\: 4\n\s+Processed reads\:\s+(?P<totalReads>\d+)\n\s+Processed bases\:\s+(?P<totalBases>\d+) bp \(.+ Mbp\)\n\s+Trimmed reads\:\s+(?P<trimmedReads>\d+) \(.+\%\)\n\s+Trimmed bases\:\s+(?P<trimmedBases>\d+) bp \(.+ Mbp\) \(.+\% of total\)\n\s+Too short reads\:\s+.+ \(.+\% of processed reads\)\n\s+Too long reads\:\s+.+ \(.+\% of processed reads\)\n\s+Total time\:\s+.+ s\n\s+Time per read\:\s+.+ ms")
		stats['qualityTrimming'][filePairId][read]            = extractData(infile=self.logPath+'/qualityTrimming.'+str(filePairId)+'.'+read+'.log.txt',           pattern='(?P<totalBasess>\d+)\tbases\n(?P<trimmedBases>\d+)\ttrimmed')

	    stats['removeEmptyReads'][filePairId] = extractData(infile=self.logPath+'/removeEmptyReads.'+str(filePairId)+'.log.txt',pattern="""Running removeEmptyReads.py:\nHeader one is empty exiting.\n(?P<totalReads>\d+) read pairs processed.\n(?P<pairsOut>\d+) read pairs to outfiles .+.\n(?P<singlets>\d+) single reads to outfile .+.\nremoveEmptyReads Exiting.""")
	    stats['bowtie2'][filePairId]          = extractData(infile=self.logPath+'/stderr.bowtie2.'+str(filePairId)+'.txt',      pattern="""(?P<totalReads>\d+) reads; of these:\n\s+(?P<pairedReads>\d+) \(\d+.\d+\%\) were paired; of these:\n\s+(?P<notPropMapedPair>\d+) \(\d+.\d+\%\) aligned concordantly 0 times\n\s+(?P<properPairs>\d+) \(\d+.\d+\%\) aligned concordantly exactly 1 time\n\s+(?P<properPairsMultiMap>\d+) \(\d+.\d+\%\) aligned concordantly >1 times\n\s+----\n\s+(?P<notPropMapedPair2>\d+) pairs aligned concordantly 0 times; of these:\n\s+(?P<discordantPairs>\d+) \(\d+.\d+\%\) aligned discordantly 1 time\n\s+----\n\s+(?P<unMappedPair>\d+) pairs aligned 0 times concordantly or discordantly; of these:\n\s+(?P<possibleSingletons>\d+) mates make up the pairs; of these:\n\s+(?P<unMappedReads>\d+) \(\d+.\d+\%\) aligned 0 times\n\s+(?P<singleSingleMap>\d+) \(\d+.\d+\%\) aligned exactly 1 time\n\s+(?P<singleMultiMap>\d+) \(\d+.\d+\%\) aligned >1 times\n(?P<overallAlignmentRate>\d+.\d+)\% overall alignment rate""")

	stats['merging'] = extractData(infile=self.logPath+'/stderr.merging.'+self.name+'.txt',pattern="Finished reading inputs.+\n.+picard.sam.MergeSamFiles done. Elapsed time",checkType='program')
	pattern = 'LIBRARY\tUNPAIRED_READS_EXAMINED\tREAD_PAIRS_EXAMINED\tUNMAPPED_READS\tUNPAIRED_READ_DUPLICATES\tREAD_PAIR_DUPLICATES\tREAD_PAIR_OPTICAL_DUPLICATES\tPERCENT_DUPLICATION\tESTIMATED_LIBRARY_SIZE\n(?P<LIBRARY>.+)\t(?P<UNPAIRED_READS_EXAMINED>\d+)\t(?P<READ_PAIRS_EXAMINED>\d+)\t(?P<UNMAPPED_READS>\d+)\t(?P<UNPAIRED_READ_DUPLICATES>\d+)\t(?P<READ_PAIR_DUPLICATES>\d+)\t(?P<READ_PAIR_OPTICAL_DUPLICATES>\d+)\t(?P<PERCENT_DUPLICATION>\d+\,\d+)\t(?P<ESTIMATED_LIBRARY_SIZE>\d+)'
	oldpattern="""LIBRARY\tUNPAIRED_READS_EXAMINED\tREAD_PAIRS_EXAMINED\tUNMAPPED_READS\tUNPAIRED_READ_DUPLICATES\tREAD_PAIR_DUPLICATES\tREAD_PAIR_OPTICAL_DUPLICATES\tPERCENT_DUPLICATION\tESTIMATED_LIBRARY_SIZE\n(?P<Library>.+)\s+(?P<unPairedReads>\d+)\s+(?P<totalReads>\d+)\s+(?P<unMapped>\d+)\s+(?P<unPairedDups>\d+)\s+(?P<pairDups>\d+)\s+(?P<opticalDups>\d+)\s+(?P<percentageDuplication>\d+\,\d+)\s+(?P<estLibSize>\d+)"""
	stats['markDuplicatesMetrix'] = extractData(infile=self.logPath+'/markDuplicatesMetrix.'+self.name+'.txt',pattern=pattern)
	stats['fixedBamFlagstat']        = extractData(infile=self.logPath+'/fixedBamFlagstat.'+self.name+'.txt',       pattern="""(?P<totalReads>\d+) \+ 0 in total \(QC-passed reads \+ QC-failed reads\)\n(?P<duplicates>\d+) \+ 0 duplicates\n(?P<mapped>\d+) \+ 0 mapped \(\d+.\d+\%:-nan\%\)\n(?P<paired>\d+) \+ 0 paired in sequencing\n(?P<read1>\d+) \+ 0 read1\n(?P<read2>\d+) \+ 0 read2\n(?P<properlyPaired>\d+) \+ 0 properly paired \(\d+.\d+\%:-nan\%\)\n(?P<bothMapped>\d+) \+ 0 with itself and mate mapped\n(?P<singletons>\d+) \+ 0 singletons \(\d+.\d+\%:-nan\%\)\n(?P<mateOnDiffChr>\d+) \+ 0 with mate mapped to a different chr\n(?P<mateOnDiffChrq5>\d+) \+ 0 with mate mapped to a different chr \(mapQ>=5\)""")
	stats['reCalibratedBamFlagstat'] = extractData(infile=self.logPath+'/reCalibratedBamFlagstat.'+self.name+'.txt',pattern="""(?P<totalReads>\d+) \+ 0 in total \(QC-passed reads \+ QC-failed reads\)\n(?P<duplicates>\d+) \+ 0 duplicates\n(?P<mapped>\d+) \+ 0 mapped \(\d+.\d+\%:-nan\%\)\n(?P<paired>\d+) \+ 0 paired in sequencing\n(?P<read1>\d+) \+ 0 read1\n(?P<read2>\d+) \+ 0 read2\n(?P<properlyPaired>\d+) \+ 0 properly paired \(\d+.\d+\%:-nan\%\)\n(?P<bothMapped>\d+) \+ 0 with itself and mate mapped\n(?P<singletons>\d+) \+ 0 singletons \(\d+.\d+\%:-nan\%\)\n(?P<mateOnDiffChr>\d+) \+ 0 with mate mapped to a different chr\n(?P<mateOnDiffChrq5>\d+) \+ 0 with mate mapped to a different chr \(mapQ>=5\)""")
	stats['unmapRemovedBamFlagstat'] = extractData(infile=self.logPath+'/unmapRemovedBamFlagstat.'+self.name+'.txt',pattern="""(?P<totalReads>\d+) \+ 0 in total \(QC-passed reads \+ QC-failed reads\)\n(?P<duplicates>\d+) \+ 0 duplicates\n(?P<mapped>\d+) \+ 0 mapped \(\d+.\d+\%:-nan\%\)\n(?P<paired>\d+) \+ 0 paired in sequencing\n(?P<read1>\d+) \+ 0 read1\n(?P<read2>\d+) \+ 0 read2\n(?P<properlyPaired>\d+) \+ 0 properly paired \(\d+.\d+\%:-nan\%\)\n(?P<bothMapped>\d+) \+ 0 with itself and mate mapped\n(?P<singletons>\d+) \+ 0 singletons \(\d+.\d+\%:-nan\%\)\n(?P<mateOnDiffChr>\d+) \+ 0 with mate mapped to a different chr\n(?P<mateOnDiffChrq5>\d+) \+ 0 with mate mapped to a different chr \(mapQ>=5\)""")
	stats['qualFilteredBamFlagstat'] = extractData(infile=self.logPath+'/qualFilteredBamFlagstat.'+self.name+'.txt',pattern="""(?P<totalReads>\d+) \+ 0 in total \(QC-passed reads \+ QC-failed reads\)\n(?P<duplicates>\d+) \+ 0 duplicates\n(?P<mapped>\d+) \+ 0 mapped \(\d+.\d+\%:-nan\%\)\n(?P<paired>\d+) \+ 0 paired in sequencing\n(?P<read1>\d+) \+ 0 read1\n(?P<read2>\d+) \+ 0 read2\n(?P<properlyPaired>\d+) \+ 0 properly paired \(\d+.\d+\%:-nan\%\)\n(?P<bothMapped>\d+) \+ 0 with itself and mate mapped\n(?P<singletons>\d+) \+ 0 singletons \(\d+.\d+\%:-nan\%\)\n(?P<mateOnDiffChr>\d+) \+ 0 with mate mapped to a different chr\n(?P<mateOnDiffChrq5>\d+) \+ 0 with mate mapped to a different chr \(mapQ>=5\)""")
	stats['noDuplicatesBamFlagstat'] = extractData(infile=self.logPath+'/noDuplicatesBamFlagstat.'+self.name+'.txt',pattern="""(?P<totalReads>\d+) \+ 0 in total \(QC-passed reads \+ QC-failed reads\)\n(?P<duplicates>\d+) \+ 0 duplicates\n(?P<mapped>\d+) \+ 0 mapped \(\d+.\d+\%:-nan\%\)\n(?P<paired>\d+) \+ 0 paired in sequencing\n(?P<read1>\d+) \+ 0 read1\n(?P<read2>\d+) \+ 0 read2\n(?P<properlyPaired>\d+) \+ 0 properly paired \(\d+.\d+\%:-nan\%\)\n(?P<bothMapped>\d+) \+ 0 with itself and mate mapped\n(?P<singletons>\d+) \+ 0 singletons \(\d+.\d+\%:-nan\%\)\n(?P<mateOnDiffChr>\d+) \+ 0 with mate mapped to a different chr\n(?P<mateOnDiffChrq5>\d+) \+ 0 with mate mapped to a different chr \(mapQ>=5\)""")
	
	#self.logPath+'/'+self.name+'.qacompute.stdout.txt'
	#self.logPath+'/'+self.name+'.qacompute.stderr.txt'
	#self.dataPath+'/'+self.name+'.qacompute.out '
	#self.logPath+'/'+self.name+'.stderr.caluclateHsmetrics.txt'
	pattern = 'READ_GROUP\n(?P<BAIT_SET>.+)\t(?P<GENOME_SIZE>\d+)\t(?P<BAIT_TERRITORY>\d+)\t(?P<TARGET_TERRITORY>\d+)\t(?P<BAIT_DESIGN_EFFICIENCY>\d+(\,\d+)?)\t(?P<TOTAL_READS>\d+)\t(?P<PF_READS>\d+)\t(?P<PF_UNIQUE_READS>\d+)\t(?P<PCT_PF_READS>\d+(\,\d+)?)\t(?P<PCT_PF_UQ_READS>\d+(\,\d+)?)\t(?P<PF_UQ_READS_ALIGNED>\d+)\t(?P<PCT_PF_UQ_READS_ALIGNED>\d+(\,\d+)?)\t(?P<PF_UQ_BASES_ALIGNED>\d+)\t(?P<ON_BAIT_BASES>\d+)\t(?P<NEAR_BAIT_BASES>\d+)\t(?P<OFF_BAIT_BASES>\d+)\t(?P<ON_TARGET_BASES>\d+)\t(?P<PCT_SELECTED_BASES>\d+(\,\d+)?)\t(?P<PCT_OFF_BAIT>\d+(\,\d+)?)\t(?P<ON_BAIT_VS_SELECTED>\d+(\,\d+)?)\t(?P<MEAN_BAIT_COVERAGE>\d+(\,\d+)?)\t(?P<MEAN_TARGET_COVERAGE>\d+(\,\d+)?)\t(?P<PCT_USABLE_BASES_ON_BAIT>\d+(\,\d+)?)\t(?P<PCT_USABLE_BASES_ON_TARGET>\d+(\,\d+)?)\t(?P<FOLD_ENRICHMENT>\d+(\,\d+)?)\t(?P<ZERO_CVG_TARGETS_PCT>\d+(\,\d+)?)\t(?P<FOLD_80_BASE_PENALTY>(\?)|(\d+?(\,\d+)?))\t(?P<PCT_TARGET_BASES_2X>(\?)|(\d+?(\,\d+)?))\t(?P<PCT_TARGET_BASES_10X>(\?)|(\d+?(\,\d+)?))\t(?P<PCT_TARGET_BASES_20X>(\?)|(\d+?(\,\d+)?))\t(?P<PCT_TARGET_BASES_30X>(\?)|(\d+?(\,\d+)?))\t(?P<PCT_TARGET_BASES_40X>(\?)|(\d+?(\,\d+)?))\t(?P<PCT_TARGET_BASES_50X>(\?)|(\d+?(\,\d+)?))\t(?P<PCT_TARGET_BASES_100X>\d+(\,\d+)?)\t(?P<HS_LIBRARY_SIZE>(\s?)|(\d+(\,\d+)?))\t(?P<HS_PENALTY_10X>\d+(\,\d+)?)\t(?P<HS_PENALTY_20X>\d+(\,\d+)?)\t(?P<HS_PENALTY_30X>\d+(\,\d+)?)\t(?P<HS_PENALTY_40X>\d+(\,\d+)?)\t(?P<HS_PENALTY_50X>\d+(\,\d+)?)\t(?P<HS_PENALTY_100X>\d+(\,\d+)?)'
	stats['hs_metrics.summary'] = extractData(infile=self.dataPath+'/'+self.name+'.hs_metrics.summary.txt',pattern=pattern)

	# make sums
	import time
	import sys
	if not list(self.getFastqs()):
	    AnalysisPipe.logfile.write('#WARNING#'+time.strftime("%Y-%m-%d:%H:%M:%S",time.localtime())+'# No fastq files found for sample: '+self.name+' continuing with next sample.\n')
	    sys.stderr.write('#WARNING#'+time.strftime("%Y-%m-%d:%H:%M:%S",time.localtime())+'# No fastq files found for sample: '+self.name+' continuing with next sample.\n')
	    self.stats = stats
	    return
	for program in ['illuminaAndNexteraTrimming','malbacWgaTrimming','qualityTrimming','rubiconWgaTrimming']:
	    try:
    		sums = {variable:0 for variable in stats[program][self.getFastqs().next()[0]]['r1'].keys()}
		for filePairId,readCount,fastq1,fastq2,sampleId in self.getFastqs():
		    for read in ['r1','r2']:
			for variable, value in stats[program][filePairId][read].iteritems(): sums[variable]+=float(value)
		stats[program]['sum']= sums
	    except AttributeError: pass
	for program in ['removeEmptyReads','bowtie2']:
	    try:
		sums = {variable:0 for variable in stats[program][self.getFastqs().next()[0]].keys()}
		for filePairId,readCount,fastq1,fastq2,sampleId in self.getFastqs():
		    for variable, value in stats[program][filePairId].iteritems(): sums[variable]+=float(value)
		stats[program]['sum']= sums
	    except AttributeError: pass
	    
	self.stats = stats

	# debug output
	#print "\n######## "+self.name+" ######## "
	#for key,value in stats.iteritems():
	#    print '    ',key
	#    try:
	#	for key2,value2 in value.iteritems():
	#	    assert type(value2) == dict
	#	    print '        ',key2,value2
	#    except:print '        ',value
	
	return 0

if __name__ == "__main__": main()