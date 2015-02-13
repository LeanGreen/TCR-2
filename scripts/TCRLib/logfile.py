#! /bin/env python

def main(): pass

def Logfile(fileName, initMessages=''):
        """ open a connection to the logfile, creates a logfile if none is present """

        # Imports
        import sys
        import time
        import os

        # check if logfile present open connection or create
        if not os.path.isfile(fileName):
            msg = '#LOGMSG#'+time.strftime("%Y-%m-%d:%H:%M:%S",time.localtime())+'# Creating the logfile "'+fileName+'".\n'
            initMessages += msg
            logFile = open(fileName,'w',1)
        else:
            logFile = open(fileName,'a',1)
            msg = '#LOGMSG#'+time.strftime("%Y-%m-%d:%H:%M:%S",time.localtime())+'# Connection to logfile '+fileName+' opened.\n'
            initMessages += msg

        logFile.write(initMessages)
        return logFile

if __name__ == "__main__": main()