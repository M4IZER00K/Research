import csv as cs
from dataclasses import fields
from re import X
import string
from tokenize import Double

proteinList = [] 
dataList = [[]]
finalData = []

with open('C:\\Users\\shash\\VS Code Programs\\Python\\ResearchWork\\metadata.csv') as csv_file:
    csv = cs.reader(csv_file, delimiter=',')
    for row in csv:
        proteinList.append(row[3])
        finalData.append([row[3]])
    finalData.pop(0)

with open('C:\\Users\\shash\\VS Code Programs\\Python\\ResearchWork\\results.csv') as csv_file:
    csv = cs.reader(csv_file, delimiter=',')
    for row in csv:
        dataList.append(row)

for dataSet in dataList[2:]:
    counter = 0
    index = 0
    for dataPoint in dataSet:
        if counter == 1:
            type = dataPoint
        if counter > 6:
            finalData[index].append((dataPoint,type))
            index += 1
        counter += 1

def compute(Protein):
    senstivity = 0.0
    specificity = 0.0
    dPos = 0.0
    dNeg = 0.0
    total = 0.0
    for dataPoint in Protein[1:]:
        if (dataPoint[1] == "Case"):
            dPos += 1
            if (float(dataPoint[0]) > 2):
                senstivity += 1
        if (dataPoint[1] == "Control"):
            dNeg += 1
            if (float(dataPoint[0]) <= 2):
                specificity += 1
    p_Y1_D1 = senstivity/dPos
    p_Y0_D0 = specificity/dNeg
    return ((p_Y1_D1,p_Y0_D0))

results = [[]]
fieldName = ["Protein","Sensitivity","Specificity"]

for protein in finalData:
    prob = compute(protein)
    results.append([protein[0],prob[0],prob[1]])
    #protein[0] -> name 
    #prob[0] -> Sensitivity
    #prob[1] -> Specificity
   
with open('resultsWeek1.csv', 'w') as csvfile:
    filewriter = cs.writer(csvfile, delimiter=',',
        quotechar='|', quoting=cs.QUOTE_MINIMAL)
    filewriter.writerow(fieldName)
    filewriter.writerows(results)
print("done")