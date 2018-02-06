##This is provcfTools

import os
import sys
import argparse
import glob
import pandas as pd
import numpy as np
import re
import datetime
from itertools import takewhile, chain

def filterAction(choices):
	#print("in action")
	class StoreDictKeyPair(argparse.Action):
		# def __init__(self, nargs=0, **kw):
		# 	super().__init__(nargs=nargs, **kw)
		def __call__(self, parser, namespace, values, option_string=None):
			#print("in call")
			my_dict={}
			for i in range(0,len(values)):
				keyVals = values[i].split("=")
				if len(keyVals)==2:
					if keyVals[0] in choices:
						vals = keyVals[1].split(",")
						my_dict[keyVals[0]] = vals
					else:
						print("ERROR: Use one of the following filter key values: "+",".join(choices))
						sys.exit(0)
				else:
					print("Error: Filter argument should be used as follows:\n--filter protein=PROT1,PROT2")
					sys.exit(0)
			setattr(namespace, self.dest, my_dict)	
	return StoreDictKeyPair
	# if not all(key in my_dict.keys for key in choices):
	# 	print("Use one of these as filter class: ["+",".join(choices)+"]")
	# 	sys.exit(1)

def filterActionCompare(choices):
	#print("in action")
	class StoreDictKeyPair(argparse.Action):
		# def __init__(self, nargs=0, **kw):
		# 	super().__init__(nargs=nargs, **kw)
		def __call__(self, parser, namespace, values, option_string=None):
			print("in call ")
			#print(values)
			my_dict={}
			print(values)
			keyVals = values
			if len(keyVals)==2:
				if keyVals[0] in choices:
					vals = keyVals[1].split(",")
					my_dict[keyVals[0]] = vals
				else:
					print("ERROR: Use one of the following Compare key values: "+",".join(choices))
					sys.exit(0)
			else:
				print("Error: Compare tool should be used as follows:\n--compare Intersection|Union|Diff vcfFile2")
				sys.exit(0)
			setattr(namespace, self.dest, my_dict)	
	return StoreDictKeyPair
	# if not all(key in my_dict.keys for key in choices):
	# 	print("Use one of these as filter class: ["+",".join(choices)+"]")
	# 	sys.exit(1)


def parseHeader(header):
	##This function parse the header list and convert this into a dictionary, key is 0th element from the list and value is 1st element.
	##Header is strore in a doctionary, all header lines will have 'string' type value except for INFO, FILTER, FORMAT
	compulsoryHDict={t[0]:t[1] for t in header if t[0]!="INFO" and t[0]!="FILTER" and t[0]!="FORMAT"}
	infoHList = [t for t in header if t[0]=="INFO"]
	filterHList = [t for t in header if t[0]=="FILTER"]
	formatHList = [t for t in header if t[0]=="FORMAT"]
	infoDF=pd.DataFrame(columns=['ID', 'Number', 'Type', 'Description', 'Source', 'Version'])
	filterDF=pd.DataFrame(columns=['ID', 'Description'])
	formatDF=pd.DataFrame(columns=['ID', 'Number', 'Type', 'Description'])
	if len(infoHList)>0:
		##generally there should be some info headers.
		temp_info=[list(chain(*(re.findall("\<ID=([^,]+),Number=([^,]+),Type=([^,]+),Description=\"([^\"]+)\"(?:,Source=\"([^\"]+)\")?(?:,Version=\"([^\"]+)\")?\>",t[1])))) for t in infoHList]
		infoDF=infoDF.append(pd.DataFrame(temp_info,columns=infoDF.columns))
	if len(filterHList)>0:
		##generally there should be some filter headers.
		temp_filter=[list(chain(*(re.findall("\<ID=([^,]+),Description=\"([^\"]+)\"\>",t[1])))) for t in filterHList]
		filterDF=filterDF.append(pd.DataFrame(temp_filter,columns=filterDF.columns))
	if len(formatHList)>0:
		##generally there should be some format headers when the provcf file is created for multiple samples.
		temp_format=[list(chain(*(re.findall("\<ID=([^,]+),Number=([^,]+),Type=([^,]+),Description=\"([^\"]+)\"\>",t[1])))) for t in formatHList]
		formatDF=formatDF.append(pd.DataFrame(temp_format,columns=formatDF.columns))
	##IDs from infoDF tells us which keys are expected in INFO column. As I propose that all the keys mentioned in INFO column
	##has to be defined in the header, provcftools should throw validation error if that is not the case.
	otherHeaders={'INFO':infoDF, 'FILTER':filterDF, 'FORMAT':formatDF}
	return [compulsoryHDict, otherHeaders]

def readVCF(filename):
	##this function read proVCF file
	with open(filename, 'r') as fobj:
	    headiter = takewhile(lambda s: s.startswith('##'), fobj)
	    header = [s.strip('\n').strip('##').split('=',maxsplit=1) for s in headiter]
	headerObj=parseHeader(header)
	vcf = pd.read_table(filename, sep='\t', keep_default_na=False, na_values=[''], header = len(header))
	###################### Validation ###############################################
	##check if all the mandatory header lines are there.
	##Mandatory headers: fileformat, source, reference, referenceSource, referenceDate, species
	mandatoryFields = ['fileformat', 'source', 'reference', 'referenceSource', 'referenceDate', 'species']
	if len(headerObj[0])>=len(mandatoryFields):
		if not all(h in headerObj[0] for h in mandatoryFields):
			##all the mandatory headers do not exist
			print("ERROR: Mandatory header missing. \nThese headers are mandatory for proVCF file: 'fileformat', 'source', 'reference', 'referenceSource', 'referenceDate', 'species'. These header names are case sensitive.")
			sys.exit(1)
	##Check in INFO column exist
	if 'INFO' in vcf.columns.values:
		if headerObj[1]['INFO'].shape[0]>0:
			#Extract key value pair from info column
			info_ext = vcf['INFO'].str.extractall("(?P<Key>[^=;]+)=(?P<Value>[^;]+)")
			##Checks if info column of every variation/row contains all the key defined in info header
			checkInfo = [headerObj[1]['INFO']['ID'][~headerObj[1]['INFO']['ID'].isin(info_ext.loc[i]['Key'])].shape[0] for i in info_ext.index.get_level_values(0).unique().tolist()]
			if max(checkInfo)>0:
				print("ERROR: INFO IDs are missing for some variations. Check variations at "+",".join([str(i+1) for i, e in enumerate(checkInfo) if e != 0]))
				sys.exit(1)
			##What happens if the keys are not in the same order???
			info_ext.index.set_levels(info_ext['Key'].unique(), level=1, inplace=True)
			info = info_ext['Value'].unstack(level=-1)
			vcf=vcf.drop('INFO',1)
			vcf=vcf.join(info)
		else:
			if len(vcf.INFO.str.split(';').tolist())>0:
				print("ERROR: INFO keys have not been defined in the header. Define all the keys in the INFO column as follows\n##INFO=<ID=ID,Number=number,Type=type,Description=\"description\",Source=\"source\",Version=\"version\"> format")
				sys.exit(1)
			else:
				print("ERROR: INFO column is empty")
				sys.exit(1)
	else:
		print("ERROR: INFO column missing. INFO column is a mandatory column of proVCF file format.")
		sys.exit(1)
	##Check FILTER column values. It should contain values that have been defined in the header, filterDF should hol all possible value for this column.
	if 'FILTER' in vcf.columns.values:
		filterValues = list(vcf['FILTER'].unique())
		if headerObj[1]['FILTER'].shape[0]>0:
			##validating filter column and header values
			if not all(f in headerObj[1]['FILTER']['ID'].tolist() for f in filterValues):
				print("ERROR: FILTER column contains undefined values. Please make sure you define all the values from FILTER column in the header.")
				sys.exit(1)
		else:
			print("ERROR: FILTER values have not been described in the header.")
			sys.exit(1)
	else:
		print("ERROR:Missing FILTER column. FILTER is a mandatory column of proVCF file format.")
		sys.exit(1)
	if 'FORMAT' in vcf.columns.values:
		formatValues = list(vcf['FORMAT'].unique())
		if headerObj[1]['FORMAT'].shape[0]>0:
			##validating format
			if not all(f in headerObj[1]['FORMAT']['ID'].tolist() for f in formatValues):
				print("ERROR: FORMAT column contains undefined values. Please make sure you define all the values from FORMAT column in the header.")
				sys.exit(1)
		else:
			print("ERROR: FORMAT IDs have not been defined in the header.")
			sys.exit(1)
	else:
		print("This file does not have format column, supposedly single sample proVcf file.")
	##Validate REF and ALT column values. These two column should only contain 21 amino acids. ALT also allows special character 'X'.
	aminoAcids=[]
	if vcf[vcf['REF'].str.contains("[^A-Z]")].shape[0]>0 or vcf[vcf['REF'].str.contains("B|J|O|U|X|Z")].shape[0]>0:
		print("ERROR: REF column contains non-amino acid letters.")
		sys.exit(1)
	if vcf[vcf['REF'].str.contains("[^A-Z]")].shape[0]>0 or vcf[vcf['ALT'].str.contains("B|J|O|U|Z")].shape[0]>0:
		print("ERROR: ALT column contains non-amino acid or 'X' letters.")
		sys.exit(1)
	##return the vcfobject and the headers
	return {'header': headerObj, 'vcf': vcf}

def writeVCF(header, vcf, fobj):
	##This function writes a proVCF file to a file name 'filename'. If it is empty, it writes to the STDOUT
	#print("writeVCF")
	
	if len(header)>0:
		##write mandatory headers:
		for k in header[0]:
			fobj.write("##"+str(k)+"="+str(header[0][k])+"\n")
		if len(header)==2:
			##Write column describing headers
			for k in header[1]:
				#print(header[1][k])
				if header[1][k].shape[0]>0:
					if k=='INFO':
						for i in range(0,header[1][k].shape[0]):
							fobj.write("##INFO=<ID="+str(header[1][k].iloc[i]['ID'])+",Number="+str(header[1][k].iloc[i]['Number'])+",Type="+str(header[1][k].iloc[i]['Type'])+",Description=\""+str(header[1][k].iloc[i]['Description'])+"\",Source=\""+str(header[1][k].iloc[i]['Source'])+"\",Version=\""+str(header[1][k].iloc[i]['Version'])+"\">\n") 
					if k=='FORMAT':
						for i in range(0,header[1][k].shape[0]):
							fobj.write("##FORMAT=<ID="+str(header[1][k].iloc[i]['ID'])+",Number="+str(header[1][k].iloc[i]['Number'])+",Type="+str(header[1][k].iloc[i]['Type'])+",Description=\""+str(header[1][k].iloc[i]['Description'])+"\">\n")
					if k=='FILTER':
						for i in range(0,header[1][k].shape[0]):
							fobj.write("##INFO=<ID="+str(header[1][k].iloc[i]['ID'])+",Description=\""+str(header[1][k].iloc[i]['Description'])+"\">\n")
	else:
		print("ERROR: vcf file is missing mandatory header lines.")
		sys.exit(1)
	#in future FORMAT should be added after checking if it exist in the header.
	fobj.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
	if vcf.shape[0]>0:
		for i in range(vcf.shape[0]):
			fobj.write(vcf.iloc[i]['#CHROM']+"\t"+str(vcf.iloc[i]['POS'])+"\t"+str(vcf.iloc[i]['ID'])+"\t"+vcf.iloc[i]['REF']+"\t"+vcf.iloc[i]['ALT']+"\t"+str(vcf.iloc[i]['QUAL'])+"\t"+vcf.iloc[i]['FILTER']+"\t")
			for j in range(header[1]['INFO']['ID'].shape[0]):
				fobj.write(str(header[1]['INFO'].iloc[j]['ID'])+"="+str(vcf.iloc[i][header[1]['INFO'].iloc[j]['ID']]))
				if j<header[1]['INFO']['ID'].shape[0]-1:
					fobj.write(";")
			fobj.write("\n")

def locationCheck(dfRow, locDF):
	test = locDF[(locDF['Protein']==dfRow['#CHROM']) & (locDF['Start'].astype(int)<=int(dfRow['POS'])) & (locDF['End'].astype(int)>=int(dfRow['POS']))]
	#print(test)
	#print(test.shape)
	if test.shape[0]>0:
		return True
	else:
		return False

def filter(vcf, value, field):
	##This function provides filtering functionality for provcftools.
	##value is considered to be a list (comma separated at the command line.)
	if field=='protein':
		##filter out all the variations reported for this protein id.
		vcf = vcf[vcf['#CHROM'].isin(value)]
	elif field == 'ref':
		##filters based on reference
		vcf = vcf[vcf['REF'].isin(value)]
	elif field == 'alt':
		##filters based on alternative sequence.
		vcf = vcf[vcf['ALT'].isin(value)]
	elif field == 'aid':
		##filters based on orf id
		if 'AID' in vcf.columns.values:
			vcf = vcf[vcf['AID'].isin(value)]
		else:
			print("There is no alternative sequence id (AID) avaliable.")
			sys.exit(1)
	elif field == 'type':
		##filters based on variation type
		if 'Type' in vcf.columns.values:
			vcf = vcf[vcf['Type'].isin(value)]
		else:
			print("There is no variation type (Type) avaliable.")
			sys.exit(1)
	elif field == 'mod':
		##modifications
		if 'Mod' in vcf.columns.values:
			vcf = vcf[vcf['Mod'].isin(value)]
		else:
			print("There is no Mod key avaliable.")
			sys.exit(1)
	elif field == 'chromosome':
		##This field appears in INFO column, if it exist, if filters out all the variations from all the proteins in this chromosome.
		if 'Chromosome' in vcf.columns.values:
			vcf = vcf[vcf['Chromosome'].isin(value)]
		else:
			print("There is no Chromosome key avaliable.")
			sys.exit(1)
	elif field == 'location':
		##location based filtering, requires protein and location information.
		##format protein_name:Start-Stop, start and stop inclusive
		locDF= pd.DataFrame(value)[0].str.extract("(?P<Protein>[^:]+):(?P<Start>\d+)-(?P<End>\d+)", expand=True)
		vcf = vcf[vcf.apply(locationCheck,args=(locDF,),axis=1)]
	elif field == 'filter':
		##Filters based on the value of filter column
		vcf = vcf[vcf['FILTER'].isin(value)]
	return vcf

def sortVcf(vcf, ascending):
	##this function will sort vcf file based on protein name/id and then varition position
	##vcf contains the vcf object, and ascending is a boolean variable to choose between ascending and descending
	return vcf.sort_values(['#CHROM','POS'], ascending=ascending)

def summarize(mHeader, vcf):
	##this function computes general summary of the provcf file.
	##Total number of variations, total proteins, total alternative sequences(if any provided)
	##Count of polymorphisms/variations by type.
	##Amino acid transition histogram for SAPs or ALTs(opt)
	##Total number of modifications and a histogram
	##Amino acid associated with modification
	print("File format: "+mHeader['fileformat'])
	print("File was generated by: "+mHeader['source'])
	print("Reference proteome: "+mHeader['reference']+" from "+mHeader['referenceSource'])
	print("Species: "+mHeader['species'])
	print("\nTotal number of variations: "+str(vcf.shape[0]))
	if 'Type' in vcf.columns.values:
		print("\tSSAPs: "+str(vcf[vcf['Type']=="SSAP"].shape[0]))
		print("\tSAPs: "+str(vcf[vcf['Type']=="SAP"].shape[0]))
		print("\tSALTs: "+str(vcf[vcf['Type']=="SALT"].shape[0]))
		print("\tALTs: "+str(vcf[vcf['Type']=="ALT"].shape[0]))
		print("\tINSs: "+str(vcf[vcf['Type']=="INS"].shape[0]))
		print("\tDELs: "+str(vcf[vcf['Type']=="DEL"].shape[0]))
	print("\nTotal Reference proteins: "+str(len(vcf['#CHROM'].unique())))
	if 'AID' in vcf.columns.values:
		print("\nTotal alternative sequences: "+str(len(vcf['AID'].unique())))
	if 'Mod' in vcf.columns.values:
		modvcf = vcf[vcf['Mod']!='.']
		mods = modvcf['Mod'].str.extract("(?P<AA>[^:]+):(?P<Position>[^:]+):(?P<Mass>.+)")
		print("Modifications reported: "+mods['Mass'].unique().to_csv())
		aggregate = {
		'Mod' : lambda x: ",".join(x)
		}
		aaMods=mods.groupby('AA').agg(aggregate)
		for i in aaMods.index:
			print(aaMods.loc[i] + " has " + ",".join(re.split(',',str(aaMods.loc[i]['Mod'])).unique()) + "modifications")
		# & ((vcf['Type']=='SAP') |(vcf['Type']=='SSAP'))



parser = argparse.ArgumentParser(description='ProVCFTools', usage='%(prog)s \n[-h\tShows available options].\n-i\tInput proVCF file[Mandatory].\n[--filter\t<FILTER_KEY>=<String>[,<String>]]\n[--sort\t[True|False]]\n[--summary]\n[--compare\t[Intersection|Union|Diff] proVCFFile2]\n[-o\tOutput proVCF File]')
parser.add_argument("-i", metavar="<FilePath>", required=True, default=sys.stdin, help="proVcf input file.")
my_dict = {}

parser.add_argument("--filter", metavar='<FILTER_KEY>=<String>[,<String...>]', nargs='+', action=filterAction(['protein','ref','alt','aid','mod','chromosome','location', 'filter']), help="\nUse any of the following filters, with one or more of them separated by comma. \
	protein=<string> : protein name \
	or ref=<string> : reference sequences \
	or alt=<string> : alternative sequences \
	or aid=<string> : alternative sequence identification \
	or mod=<string> : modification name or mass shift. If mass shift is used, put in a quote \
	or chromosome=<string> : chromosome number\
	or location=<string> : Each location has protein_id:start_location-end_location format \
	or type=<string> : variation type, e.g. SAP,INS,DEL \
	or filter=<string> : Filters based on the filters applied to the data.")
parser.add_argument("--sort", metavar='[True|False]', nargs=1, choices=['True','False'],help="Sort the input file. Use --sort True for ascending order and --sort False for descending order.")
parser.add_argument("--summary", action='store_false', help="A summary of the input file.")
parser.add_argument("--compare", metavar=('[Intersection|Union|Diff]', '[<FILE_PATH>,<FILE_PATH>]'), nargs=2, action=filterActionCompare(['Intersection','Union','Diff']), help="\nCompare input files and return intersection or union of variations of these input files. Use Intersection or Union or Difference (variations only appearing in one file) and provide list of files you want to compare.")
parser.add_argument("-o", metavar="<FilePath>", default=sys.stdout, help="proVcf output file.")
args = parser.parse_args()
print("Input File:" + args.i)
#print("Summary File:" + args.summary)
vcfDict=readVCF(args.i)
vcf=vcfDict['vcf']
if args.o!=sys.stdout:
	fobj=open(args.o, 'w')
else:
	fobj=sys.stdout
#print("Shape of initial vcf obj:"+str(vcfDict['vcf'].shape[0]))
if args.summary==True:
	#print("Summary File:" + str(args.summary))
	summarize(vcfDict['header'][0], vcfDict['vcf'])
else:
	print("Summary not requested")
if args.filter!=None:
	for i in args.filter:
		vcf = filter(vcfDict['vcf'], args.filter[i], i)
		print("Shape of vcf obj after applying each filter:"+str(vcf.shape[0]))
if args.sort!=None:
	if args.sort[0]=='True':
		vcf=sortVcf(vcf, True)
		#print("Shape of vcf obj:"+str(vcf.shape[0]))
		#writeVCF(vcfDict['header'],vcf,fobj)
	elif args.sort[0]=='False':
		vcf=sortVcf(vcf, False)
		print("Shape of vcf obj:"+str(vcf.shape[0]))
writeVCF(vcfDict['header'],vcf,fobj)