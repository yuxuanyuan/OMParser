#!/usr/bin/env python3

# This script aims to parse the OM alignments from OMBlast and get the coverage 
# ScriptName: omCoverage.py
# Created date: 22/04/2020
# Last modified: 22/04/2020
# Notes: v1.0 is the initial version
# Copyright: Copyright (c) 2020 Yuxuan Yuan (yuxuan.yuan@outlook.com)

#======================= Modules ======================
import sys
import os
import logging
import argparse
import numpy as np
import pandas as pd

#======================== Function ==========================
DEBUG="" #change it when debugging
logFormat = "%(asctime)s [%(levelname)s] %(message)s"
level = "DEBUG" if DEBUG != "" else "INFO"
logging.basicConfig( stream=sys.stderr, level=level, format=logFormat )

cmapHeader = "# CMAP File Version:\t0.1\n\
# Label Channels:\t1\n\
# Nickase Recognition Site 1:\tunkbown\n\
# Number of Consensus Maps:\tunknown\n\
#h CMapID\tContigLength\tNumSites\tSiteID\tLabelChannel\tPosition\tStdDev\tCoverage\tOccurrence\n\
#f int\tfloat\tint\tint\tint\tfloat\tfloat\tint\tint\
"

class omCov(object):
    """check the OM coverage of chromosome specific regions"""
    
    def __init__(self, oma, refcmap, chr, start, end, prefix, outDir):
        self.oma = oma
        self.refcmap = refcmap
        self.chr = chr
        self.start = int(start)
        self.end = int(end)
        self.prefix = prefix
        self.outDir = outDir

    def checkDir(self):
        """check the output directory"""
        logging.info("Checking folder: '%s'" % self.outDir)
        dirname = os.path.abspath(self.outDir)
        if not os.path.isdir(dirname):
            logging.error("Folder: '%s' does not exit. Please check!" % self.outDir)
            sys.exit(-1)
        if not os.access(dirname, os.W_OK):
            logging.error("Folder folder: '%s' is not writable. Please check!" % self.outDir)
            sys.exit(-1) 

    def checkFile(self, Filename):
        """check the input file"""
        logging.info("Checking file: '%s'" % Filename) 
        filename =  os.path.abspath(Filename)
        if not os.path.isfile(filename):
            logging.error("File: '%s' does not exit. Please check!" % Filename)
            sys.exit(-1)
        if not os.access(filename, os.R_OK):
            logging.error("File: '%s' is not readable. Please check!" % Filename)
            sys.exit(-1) 

    def parseAlignment(self):
        """parse the alignment"""
        logging.info("Starting to parse the alignemnt ...")
        refCoords = dict()
        qryCoords = dict()
        parsed = np.empty(shape = [0, 8])

        logging.info("Loading the refcmap file ...")
        ## parse refcmap
        with open (self.refcmap, 'r') as ref:
            for line in ref:
                line = line.strip()
                if line[0] != "#":
                    items = line.split("\t")
                    refID = items[0]
                    #seqLen = int(float(items[1]))
                    nSites = int(items[2])
                    siteID = int(items[3])
                    siteCoord = int(float(items[5]))
                    if siteID <= nSites:
                        refCoords['%s-%s'%(refID, siteID)] = siteCoord
        
        logging.info("Parsing the OMA file ...")
        ##parse the OM file
        with open (self.oma, 'r') as OMA:
            for line in OMA:
                sum = 0
                line = line.strip()
                if line[0] != "#":
                    items = line.split("\t")
                    molID = items[0]
                    #nMolSeq = int(items[1])
                    molSegs = items[2]
                    coords = molSegs.split(";")
                    S = 0.0
                    for i in range(len(coords)):
                        S += int(coords[i])
                    S += len(coords) - 1
                    for i in range(len(coords)):                                
                        sum = 0.0
                        for j in range(i+1):
                            sum += int(coords[j])
                        qryCoords['%s-%s'%(molID, i+1)] = sum+i+1

                    try:
                        refID = items[3]
                        strand = items[4]
                        score = float(items[5])
                        #confidence = float(items[6])
                        refSegStartID = int(items[7])
                        refSegEndID = int(items[8])
                        molSegStartID = int(items[9])
                        molSegEndID = int(items[10])
                        #refStartCoord = int(items[11])
                        #refEndCoord = int(items[12])
                        cigar = items[13]

                        mylist = []
                        start = 0
                        i = 0 
                        while  i <len(cigar):
                            try: 
                                int(cigar[i])
                            except ValueError:
                                if cigar[i] == 'M':
                                    match = int(cigar[start:i])
                                    mylist.append((match, 'M'))
                                if cigar[i] == 'D':
                                    deletion = int(cigar[start:i])
                                    mylist.append((deletion, 'D'))
                                if cigar[i] == 'I':
                                    insertion = int(cigar[start:i])
                                    mylist.append((insertion, 'I'))
                                start = i + 1
                            i += 1

                        if strand == "+":
                            for i in range(len(mylist)):
                                if mylist[i][1] == 'M':
                                    j = 0
                                    while j < int(mylist[i][0]):
                                        try:
                                            refPos = refCoords['%s-%s'%(refID, refSegStartID+j)]
                                            molPos = qryCoords['%s-%s'%(molID, molSegStartID+j)]
                                            if molSegStartID +j <= molSegEndID+1:
                                                parsed = np.append(parsed,[[refID, refSegStartID+j, refPos, molID, molSegStartID+j, molPos, strand, score]], axis = 0)
                                        except KeyError:
                                            pass
                                        j += 1
                                    refSegStartID += j
                                    molSegStartID += j

                                if mylist[i][1] == 'D':
                                    refSegStartID +=  int(mylist[i][0])
                                if mylist[i][1] == 'I':
                                    molSegStartID += int(mylist[i][0])                
                        if strand == "-":
                            molSegStartID += 1
                            i = 0
                            while i < len(mylist):
                                if mylist[i][1] == 'M':
                                    j = 0
                                    while j < int(mylist[i][0]):
                                        try:
                                            refPos = refCoords['%s-%s'%(refID, refSegStartID+j)]
                                            molPos = qryCoords['%s-%s'%(molID, molSegStartID-j)]
                                            if refSegStartID <= refSegEndID+1:
                                                parsed = np.append(parsed, [[refID, refSegStartID+j, refPos, molID, molSegStartID-j, molPos, strand, score]], axis = 0) 
                                        except KeyError:
                                            pass
                                        j += 1
                                    refSegStartID += j
                                    molSegStartID -= j     
                                                                        
                                if mylist[i][1] == 'D':
                                    refSegStartID += int(mylist[i][0])
                                if mylist[i][1] == 'I':
                                    molSegStartID -=  int(mylist[i][0]) 
                                i += 1
                    except IndexError:
                        pass
        self.parsed = pd.DataFrame(parsed, columns=['RefID', 'RefSiteID', 'RefSitePos', 'QryID', 'QrySiteID', 'QrySitePos', 'Strand', 'Score'])
        #self.parsed.to_csv("%s/%s.parsed.txt" % (self.outDir, self.prefix), sep = '\t', index = False, header = True)
        logging.info("Completing the parsing ...")

    def getCov(self):
        logging.info("Writing the coverage information to '%s.cov.txt'"% self.prefix)
        with open('%s/%s.cov.txt' % (self.outDir, self.prefix), 'w') as cov:
            cov.write('RefID\tRefSiteID\tRefSitePos\tCoverage\n')
            tmp = self.parsed[self.parsed['RefID'] == self.chr].reset_index(drop=True)
            if len(tmp) == 0:
                logging.error("Oops! It seem the input chrID: '%s' is not correct! Please check!" % self.chr)
                sys.exit(-1)
            if self.end < self.start:
                logging.error("Oops! It seem the end postion is smaller the start position! Please check!")
                sys.exit(-1)                
            sortTmp = tmp.sort_values(['RefSitePos']).reset_index(drop=True)
            sortTmp = sortTmp.astype({'RefSitePos': 'int32'})
            target = sortTmp[(sortTmp['RefSitePos']<=self.end) & (sortTmp['RefSitePos']>=self.start)]
            for i in target['RefSiteID'].unique():
                sub = target[target['RefSiteID'] == i].reset_index(drop=True)
                coverage = len(sub)
                pos = sub['RefSitePos'][0]
                cov.write('%s\t%s\t%s\t%s\n' % (self.chr, i, pos, coverage))


if __name__=="__main__":
    """
    Parse the function and usage.
    """
    parser = argparse.ArgumentParser(description='Parse the OM alignments from OMBlast')
    parser.add_argument('-v', '--version', action='version', version='1.0')
    parser.add_argument('-O', dest='oma', help='OMA file from OMBlast', type = str)
    parser.add_argument('-r', dest='refCMAP', help='reference cmap file', type = str) 
    parser.add_argument('-c', dest='chr', help='the chromosome ID', type = str)   
    parser.add_argument('-s', dest='start', help='start position', default = 0, type=int)
    parser.add_argument('-e', dest='end', help='end position', default = 0, type=int)
    parser.add_argument('-p', dest='prefix', help='prefix of the output', type=str)
    parser.add_argument('-o', dest='outDir', help='output directory', type=str)
    args = parser.parse_args()
    if None not in [args.oma, args.refCMAP, args.chr, args.start, args.end, args.prefix, args.outDir]:            
        run = omCov(args.oma, args.refCMAP, args.chr, args.start, args.end, args.prefix, args.outDir)
        run.checkDir()
        run.checkFile(args.oma)
        run.checkFile(args.refCMAP)
        run.parseAlignment()
        run.getCov()
    else:
        print
        parser.print_help()
        print
        sys.exit(-1)
