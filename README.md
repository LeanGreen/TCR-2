# TCR
A pipe for extracting TCRs from RNAseq data.  

To start a new analysis run the following commands

1. Optionally start with:  
 `scripts/TCRcaller init <pathToAnalysis>`  
to create folders and database structures logfile etc.

2. Then run the add sample command to add a sample id/name
 `scripts/TCRcaller addSample <pathToAnalysis> <sampleName>`  

3. When a sample has been added it's time to add some data to that sample either by adding a pair of fastq files by the following command:
`scripts/TCRcaller addFastq  <pathToAnalysis> <sampleName> <read1file> <read2file>`  
or by adding a single fastq file:  
`scripts/TCRcaller addFastq  <pathToAnalysis> <sampleName> <singleReadsFile>` 

4. If its the first time the analysis is run on the system you will need to download the tools needed this can be done by adding the term `getTools` anywhere on the commandline eg:  
`scripts/TCRcaller getTools <pathToAnalysis> <sampleName>`

5. Additionally a bowtie2 reference with TCR sequences is needed for the alignment to work. The software automatically uses the latest created/changed reference file matching the following pattern `references/*.fa`. Either you create your own reference or add the term `getReferences` anywhere on the commandline to download the latest sequences from NCBI and build a bowtie2 index. For example like this:  
`scripts/TCRcaller getTools <pathToAnalysis> <sampleName> getReferences`  
to both download tools and create a bowtie2 reference in one step.

6. Finally to tell the software to do some work run the following command:  
`scripts/TCRcaller findTCR <pathToAnalysis>`

Have fun!!  
and tell me what you want next ;)  
//EB