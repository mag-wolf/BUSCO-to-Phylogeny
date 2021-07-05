#!/usr/bin/python3
import os, sys, os.path

path = sys.argv[1]  #open stuff


# Open a file
dirs = os.listdir( path )

listS = []
listA = []
listO = []
keepS = {}
counter = 0

for content in os.listdir(path):
	if content.startswith('run_'):
		counter = counter + 1

exludecount = counter // 4
finalcout = counter - exludecount
#print(finalcout)


for content in os.listdir(path):     
	if content.startswith('run_'):
		pathscos=path + "/" + content + "/" + "busco_sequences/single_copy_busco_sequences/"
		for file in os.listdir(pathscos):
			if file.endswith(".faa"):
				listA.append(file)
				if file not in keepS:
					keepS[file] = [0]
					listS.append(file)

for itemS in listS:
	countS=listA.count(itemS)
#	print(itemS)
#	print(countS)
	if countS >= finalcout:
		listO.append(itemS)

for itemO in listO:
	for content in os.listdir(path):
		if content.startswith('run_'):
			contentname=content.replace("run_", "")
			pathscos2=path + "/" + content + "/" + "busco_sequences/single_copy_busco_sequences/"
			itemOsplit,rest = itemO.split(".",1)
			filenameS = os.path.join(pathscos2, itemO)
#			#print(filenameS)
			if os.path.exists(filenameS):
#				#print (filenameS, " exists")			
				#filenameS = os.path.join(pathscos2, itemO)
				filehandleS = open(filenameS)
				input=filehandleS.read()
				input2=input.split("\n",1)[1];
				input3=">" + contentname + " " + "\n" + input2
				filenameN = itemOsplit + "_complete.fasta"
				filehandleN = open(filenameN,"a")
				filehandleN.write(input3)
				filehandleS.close()
			#	filehandleN.close()

