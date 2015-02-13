#! /bin/env python

def main(): pass

def test():
    print "tjoheh"

class Database(object):

    def __init__(self, dbPath, TCRcaller, logfile=''):
	self.path = dbPath
	self.logfile=logfile
	self.TCRcaller=TCRcaller
	if type(self.logfile) == str():
	    import sys
	    sys.stderr.write('LOGFILE ERROR! IN DTATBASE')
	    sys.exit()

    def getConnection(self,):
	#
	# Import useful stuff
	#
	import sqlite3
	import sys

	#
	# Create database and set
	#
	try: self.conn = sqlite3.connect(self.path)
	except sqlite3.OperationalError:
	    print 'ERROR: Trouble with the database, plase check your commandline.'
	    sys.exit()
	self.c = self.conn.cursor()

    def commitAndClose(self,):
	#
	# commit changes and close connection
	#
	self.conn.commit()
	self.conn.close()

    def create(self,):
	""" creates the database holding all information used in the analysis """
	
	self.getConnection()
	
	#
	# Create tables
	#
	self.c.execute('''CREATE TABLE runs (startTime,command,commandLine,finishedSuccessfully,masterPid)''')
	self.c.execute('''CREATE TABLE fastqs (filePairId,fastq1,fastq2,readCount,addedToReadsTable,minReadLength,sampleId,PRIMARY KEY (filePairId))''');
	self.c.execute('''CREATE TABLE settings (variableName,defaultValue,value,setTime,PRIMARY KEY (variableName))''')
	self.c.execute('''CREATE TABLE results (resultName,defaultValue,value,setTime,PRIMARY KEY (resultName))''')
	self.c.execute('''CREATE TABLE samples (sampleId,sampleName,refType,PRIMARY KEY (sampleId))''')
	
	self.commitAndClose()

    def addToRunsTable(self, startTime, command, commandLine, finishedSuccessfully, masterPid):
	
	self.getConnection()
	
	#
	# check if pid already in database
	#
	t = (masterPid,)
	data = self.c.execute('SELECT masterPid, startTime FROM runs WHERE masterPid=?',t).fetchall()        
	if data:
	    for tmp1,tmp2 in data:

	#
	# if pid and startTime matches update the "finishedSuccessfully" entry
	#
		if tmp1 == masterPid and tmp2 == startTime:
		    values = (startTime, command, commandLine, finishedSuccessfully, masterPid)
		    self.c.execute('UPDATE runs SET finishedSuccessfully=? WHERE masterPid=? AND startTime=?', (finishedSuccessfully,masterPid,startTime))
	
	#
	# if not in the database add a new row
	#
	else:
	    values = (startTime, command, commandLine, finishedSuccessfully, masterPid)
	    self.c.execute('INSERT INTO runs VALUES (?,?,?,?,?)', values)
	
	self.commitAndClose()
	
	return 0

    def addSample(self, newSampleName,newSampleRefType=None):
	
	#
	# Imports
	#
	import sys
	import time
	import sample
	
	#
	# open connection to database
	#
	self.getConnection()
	
	sampleNames = []
	sampleIds = []
	
	#
	# check if any of the fastqs already in database
	#
	data = self.c.execute('SELECT sampleId,sampleName,refType FROM samples').fetchall()
	if data:
	    for (sampleId,sampleName,sampleRefType) in data:
		#sampleName = sampleName[0]
		sampleNames.append(sampleName)
		sampleIds.append(sampleId)
	    if newSampleName in sampleNames:
		msg = '#ERROR_MSG#'+time.strftime("%Y-%m-%d:%H:%M:%S",time.localtime())+'# SampleName must be uniq, there is already a sample with name '+newSampleName+' , exiting.\n'
		self.logfile.write(msg)
		sys.stderr.write(msg)
		sys.exit(1)

	
	if sampleIds:  sampleId = max(sampleIds)+1
	else:          sampleId = 0 
	self.logfile.write('#LOGMSG#'+time.strftime("%Y-%m-%d:%H:%M:%S",time.localtime())+'# Adding sample '+newSampleName+' to database with id '+str(sampleId)+'.\n')
	values = (sampleId,newSampleName,newSampleRefType)
	self.c.execute('INSERT INTO samples VALUES (?,?,?)', values)
	
	sample = sample.Sample(self.TCRcaller, sampleName=newSampleName, sampleId=sampleId,refType=newSampleRefType)
	sample.createDirs()
	
	self.commitAndClose()
	
	return 0

    def getSamples(self):
	#
	# Imports
	#
	import sys
	import time
	
	#
	# open connection to database
	#
	self.getConnection()
	
	refSamples = []
	samples = []
	
	data = self.c.execute('SELECT sampleId,sampleName,refType FROM samples').fetchall()
	if data:
	    for (sampleId,sampleName,sampleRefType) in data:
		if sampleRefType: refSamples.append( Sample(sampleName=sampleName,sampleId=int(sampleId),refType=sampleRefType) )
		else:                samples.append( Sample(sampleName=sampleName,sampleId=int(sampleId),refType=None) )
	
	self.commitAndClose()
	
	return refSamples+samples

    def updateFastqReadCount(self,sample):

	self.getConnection()

	#
	# check if any of the fastqs already in database
	#
	data = self.c.execute('SELECT filePairId,fastq1,fastq2,readCount,addedToReadsTable,minReadLength,sampleId FROM fastqs').fetchall()
	if data:
	    for filePair in data:
		if int(sample.id) == int(filePair[-1]):
		    filePairId,fastq1,fastq2,readCount,addedToReadsTable,minReadLength,sampleId = filePair
		    tmp = extractData(infile=sample.logPath+'/rubiconWgaTrimming.'+str(filePairId)+'.r1.log.txt',        pattern="Running wgaAdapterTrimmer.py\nProcessed a total of\t(?P<totalReads>\d+)\treads. \(.+\)\nProcessed a total of\t(?P<totalBases>\d+)\tbases \(.+\).\ntrimmed a total of\t(?P<trimmedBases>\d+)\tbases in the start of reads \(.+\).\nwgaAdapterTrimmer.py done exiting ...\n?")
		    if type(tmp) != str:newreadcount = int(tmp['totalReads'])
		    else: newreadcount = 'Unknown'
		    self.c.execute('UPDATE fastqs SET readCount=? WHERE filePairId=?', (newreadcount,filePairId))

	self.commitAndClose()

    def addFastq(self, sampleNameOrId, fastq1, fastq2):

	#
	# Imports
	#
	import sys
	import os
	import time
	
	fastq1 = os.path.abspath(fastq1)
	fastq2 = os.path.abspath(fastq2)
	
	samples = TCRcaller.database.getSamples()
	samplesbyName = {}
	samplesbyId = {}
	for sample in samples:
	    samplesbyId[sample.id]=sample
	    samplesbyName[sample.name]=sample
	sampleName = None
	sampleId = None
	try:
	    if sampleNameOrId.isdigit() and (int(sampleNameOrId) in [sample.id for sample in samples]):
		sampleId = int(sampleNameOrId);
		sampleName = str(samplesbyId[int(sampleId)].name)
	except ValueError: pass
	if type(sampleId) == int and type(sampleName) == str: pass
	elif   sampleNameOrId  in samplesbyName.keys():
	    sampleName = sampleNameOrId;
	    sampleId = samplesbyName[sampleName].id
	else:
	    msg = '#ERROR_MSG#'+time.strftime("%Y-%m-%d:%H:%M:%S",time.localtime())+'# SampleName (or id) must be registered in the database, there is no sample with name or id '+str(sampleNameOrId)+' ('+str(type(sampleNameOrId))+') , exiting.\n'
	    self.logfile.write(msg)
	    sys.stderr.write(msg)
	    sys.exit(1)

	#
	# open connection to database
	#
	self.getConnection()
	
	filePairId = None
	filePairIds = []
	
	#
	# check if any of the fastqs already in database
	#
	data = self.c.execute('SELECT filePairId,fastq1,fastq2 FROM fastqs').fetchall()
	if data:
	    for filePair in data:
		filePairId = int(filePair[0])
		filePairIds.append(filePairId)
		for fastq in [fastq1, fastq2]:
		    if fastq in filePair:
			message = 'ERROR: '+fastq+' already in the database.\nExiting after error.'
			print message
			self.logfile.write(message+'\n')
			sys.exit(1)
	#
	# if not in the database add a new row
	#
	self.logfile.write('Getting readcount for file'+fastq1+' ... \n')
	readCount = 'Unknown'#bufcount(fastq1)/4 #one read is four lines
	self.logfile.write('...done. The file has '+str(readCount)+' reads.\n')
	addedToReadsTable = False#SEAseqPipeLine.startTimeStr
	minReadLength = 'NA'

	if filePairIds: filePairId = max(filePairIds)+1
	else: filePairId = 0
	values = (filePairId,fastq1,fastq2,readCount,addedToReadsTable,minReadLength,sampleId)
	self.c.execute('INSERT INTO fastqs VALUES (?,?,?,?,?,?,?)', values)
	
	self.commitAndClose()
	
	return 0

    def getFastqs(self,):
	#
	# Imports
	#
	import sys
	
	#
	# open connection to database
	#
	self.getConnection()
		
	#
	# get att data in fastqs table
	#
	filePairs = self.c.execute('SELECT filePairId,fastq1,fastq2,readCount,addedToReadsTable,minReadLength,sampleId FROM fastqs').fetchall()
	
	self.commitAndClose()
	
	#return [[readCount,fastq1,fastq2] if (not addedToReadsTable) else None for filePairId,fastq1,fastq2,readCount,addedToReadsTable,minReadLength in filePairs]
	return [[filePairId,readCount,fastq1,fastq2,sampleId] for filePairId,fastq1,fastq2,readCount,addedToReadsTable,minReadLength,sampleId in filePairs]

    def getRuns(self, runTypes):

	self.getConnection()

	runsInfo = []
	data = self.c.execute('SELECT * FROM runs').fetchall()
	for startTime, command, commandLine, finishedSuccessfully, masterPid in data:
	    if command in runTypes: runsInfo.append([startTime, command, commandLine, finishedSuccessfully, masterPid])

	self.commitAndClose()

	return runsInfo

if __name__ == "__main__": main()