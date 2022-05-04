import csv as cs
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
    total = 0.0
    for dataPoint in Protein[1:]:
        if (float(dataPoint[0]) > 2 or dataPoint[1] == "Case"):
            senstivity += 1
        if (float(dataPoint[0]) <= 2 or dataPoint[1] == "Control"):
            specificity += 1
        total += 1
    p_Y1_D1 = senstivity/total
    p_Y0_D0 = specificity/total
    return ((p_Y1_D1,p_Y0_D0))

for protein in finalData:
    prob = compute(protein)
    print(protein[0])
    print("Sensitivity = ", prob[0])
    print("Specificity = ", prob[1])
    print("===========================================")
    print()