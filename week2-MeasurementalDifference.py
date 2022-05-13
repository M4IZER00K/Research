#PARAMATIZE THE FUNCTION TO COMBINE RELATED SAMPLES OF BLOOD MARKERS 

from cgi import test
from contextlib import nullcontext
import csv as cs
from dataclasses import fields
from operator import indexOf
from queue import PriorityQueue
from re import X, sub
import string
from this import d
from tokenize import Double

week1Data = [[]]
dataList = [[]]
finalData = [[]]
S2MetricData = [[]]
subsetData = [[]]
subsetD_r_ = [[]]
subsetR_r = [[]]
ProteinData = {"test": []}

conjunctiveData = {"test": float}
conjunctiveData.pop("test")

disjunctiveData = {"test": float}
disjunctiveData.pop("test")

finalReturnData = []

averageS2 = 0.0
with open('C:\\Users\\shash\\VS Code Programs\\Python\\ResearchWork\\resultsWeek1.csv') as csv_file:
    csv = cs.reader(csv_file, delimiter=',')
    week1Data.pop(0)
    for row in csv:
        if len(row) != 0:
            week1Data.append(row)
            finalData.append([row[0]])
    week1Data.pop(0)
    finalData.pop(0)
    finalData.pop(0)

with open('C:\\Users\\shash\\VS Code Programs\\Python\\ResearchWork\\results.csv') as csv_file:
    csv = cs.reader(csv_file, delimiter=',')
    for row in csv:
        dataList.append(row)

#correlates all the data with its Protein
for dataSet in dataList[2:]:
    counter = 0
    index = 0
    for dataPoint in dataSet:
        if counter == 0:
            ig = dataPoint
        elif counter == 1:
            type = dataPoint
        elif counter == 2:
            stage = dataPoint
        elif counter == 4:
            age = dataPoint
        if counter > 6:
            finalData[index].append((dataPoint,ig,type,stage,age))
            index += 1
        counter += 1

#creates Dict for Protein : data 
for Data in finalData:
    ProteinData[Data[0]] = Data[1:]
S2MetricData.pop(0)
avg = 0.0
count = 0.0
for Protein in week1Data:
    s2 = float(Protein[1])*float(Protein[2])
    S2MetricData.append([Protein[0],s2])
    avg += s2
    count += 1
averageS2 = avg / count

#create subset of markers that have s2 above average
S2MetricProteinKey = {"test": float}
S2MetricProteinKey.pop("test")
for Protein in S2MetricData:
    if (Protein[1] > averageS2):
        subsetData.append(Protein[0])
        S2MetricProteinKey[Protein[0]] = Protein[1]
subsetData.pop(0)
print(len(subsetData)) 

#D_r_ has each protein with the case and the c value, while D_r_Cases has each protein with only the case
D_r_ = {"test": []}
D_r_.pop("test")
D_r_Cases = {"test": []}
D_r_Cases.pop("test")
for Protein in subsetData:
    data = ProteinData[Protein]
    tempData = []
    tempDataCases = []
    for point in data:
        if (float(point[0]) > 2 ):
            tempData.append(point)
            tempDataCases.append(point[1:])
    D_r_[Protein] = tempData
    D_r_Cases[Protein] = tempDataCases
print(*D_r_Cases['HpCD00780232'])

#given input of two Proteins, checks for similiarity of cases and depending on minimumSamples (variable to change) -> returns True if enough similar cases and False otherwise
def Similiarity(Protein1, Protein2):
    dataSet1 = D_r_Cases[Protein1]
    dataSet2 = D_r_Cases[Protein2]
    sharedSamples = 0
    minimumSamples = 1
    for data in dataSet1:

        if data in dataSet2:
            sharedSamples += 1
    if sharedSamples >= minimumSamples:
        return True
    else:
        return False

#computes the S2 Metric for each case
def compute(data):
    senstivity = 0.0
    specificity = 0.0
    dPos = 0.0
    dNeg = 0.0
    for dataPoint in data:
        if (dataPoint[1] == "Case"):
            dPos += 1
            if (float(dataPoint[0]) == 1.0):
                senstivity += 1
        if (dataPoint[1] == "Control"): 
            dNeg += 1
            if (float(dataPoint[0]) == 0.0):
                specificity += 1
    p_Y1_D1 = senstivity/dPos
    p_Y0_D0 = specificity/dNeg
    return (p_Y1_D1*p_Y0_D0)

#only finds the values for the first 500 out of ~1014 Proteins whose S2 Metric > avergage
def conjunctive():
    R_r_ = {"test": []}
    notR_r_ = {"test": []}
    R_r_.pop("test")
    notR_r_.pop("test")
    for Protein in subsetData[:500]:
        tempProteinCorrelatedData = []
        notTempProteinCorrelatedData = []
        for SecondProtein in subsetData:
            if (Similiarity(Protein, SecondProtein) and Protein != SecondProtein):
                dataset1 = ProteinData[Protein]
                dataset2 = ProteinData[SecondProtein]
                tempSampleData = []
                notTempSampleData = []
                for index in range(len(dataset1)):
                    test = float(dataset1[index][0])
                    if (float(dataset1[index][0]) > 2 and float(dataset2[index][0]) > 2):
                        tempSampleData.append((1,dataset1[index][2]))
                    else:
                        tempSampleData.append((0, dataset1[index][2]))
                    if (float(dataset1[index][0]) <= 2 and float(dataset2[index][0]) > 2):
                        notTempSampleData.append((1,dataset1[index][2]))
                    else:
                        notTempSampleData.append((0,dataset1[index][2]))
                tempProteinCorrelatedData.append((SecondProtein,tempSampleData))
                notTempProteinCorrelatedData.append((SecondProtein,notTempSampleData))
        if (len(tempProteinCorrelatedData) > 0):
            R_r_[Protein] = tempProteinCorrelatedData
        if (len(notTempProteinCorrelatedData) > 0):
            notR_r_[Protein] = notTempProteinCorrelatedData
    #finds Eavg for the Protein
    for Protein in R_r_:
        eavg = 0
        for index in range(len(R_r_[Protein])):
            dataP_r_r = R_r_[Protein][index]
            datanotP_r_r = notR_r_[Protein][index]
            Pr_r = compute(dataP_r_r[1])
            Pnotr_r = compute(datanotP_r_r[1])
            eavg += Pr_r - Pnotr_r
        eavg = eavg/len(R_r_[Protein])
        conjunctiveData[Protein] = eavg

#only finds the values for the first 500 out of ~1014 Proteins whose S2 Metric > avergage
def disjunctive():
    R_r_ = {"test": []}
    notR_r_ = {"test": []}
    R_r_.pop("test")
    notR_r_.pop("test")
    for Protein in subsetData[:500]:
        tempProteinCorrelatedData = []
        notTempProteinCorrelatedData = []
        for SecondProtein in subsetData:
            if (Similiarity(Protein, SecondProtein) and Protein != SecondProtein):
                dataset1 = ProteinData[Protein]
                dataset2 = ProteinData[SecondProtein]
                print(Protein, "  ", SecondProtein)
                tempSampleData = []
                notTempSampleData = []
                for index in range(len(dataset1)):
                    test = float(dataset1[index][0])
                    if (float(dataset1[index][0]) > 2 or float(dataset2[index][0]) > 2):
                        tempSampleData.append((1,dataset1[index][2]))
                    else:
                        tempSampleData.append((0, dataset1[index][2]))
                    if (float(dataset1[index][0]) <= 2 or float(dataset2[index][0]) > 2):
                        notTempSampleData.append((1,dataset1[index][2]))
                    else:
                        notTempSampleData.append((0,dataset1[index][2]))
                tempProteinCorrelatedData.append((SecondProtein,tempSampleData))
                notTempProteinCorrelatedData.append((SecondProtein,notTempSampleData))
        if (len(tempProteinCorrelatedData) > 0):
            R_r_[Protein] = tempProteinCorrelatedData
        if (len(notTempProteinCorrelatedData) > 0):
            notR_r_[Protein] = notTempProteinCorrelatedData
    for Protein in R_r_:
        eavg = 0
        for index in range(len(R_r_[Protein])):
            dataP_r_r = R_r_[Protein][index]
            datanotP_r_r = notR_r_[Protein][index]
            Pr_r = compute(dataP_r_r[1])
            Pnotr_r = compute(datanotP_r_r[1])
            eavg += Pr_r - Pnotr_r
        eavg = eavg/len(R_r_[Protein])
        disjunctiveData[Protein] = eavg
    return disjunctiveData

def finalReturn():
    S2Data = S2MetricProteinKey

    # "n/a" means that the specific Protein didnt have a value for conjunctive or disjunctive because it didn't share any common blood markers
    for Protein in S2Data:
        conjunctiveValue = ""
        disjunctiveValue = ""
        if Protein in conjunctiveData:
            print("run")
            print(Protein)
            print(conjunctiveData[Protein])
            conjunctiveValue = conjunctiveData[Protein]
            print(conjunctiveValue)
        else:
            conjunctiveValue = "n/a"
        if Protein in disjunctiveData:
            disjunctiveValue = disjunctiveData[Protein]
        else:
            disjunctiveValue = "n/a"
        finalReturnData.append((Protein,S2MetricProteinKey[Protein],conjunctiveValue,disjunctiveValue))
        #finalReturnData = ("Protein Name",S2 Metric, Conjuctive Eavg, Disjunctive Eavg)

def writeCSV():
    fieldName = ["Protein","S2 Metric", "Conjunctive Eavg","Disjunctive Eavg"]

    with open('resultsWeek2(S2_Metric_Conjuctive_Disjunctive_Eavg).csv', 'w') as csvfile:
        filewriter = cs.writer(csvfile, delimiter=',',
            quotechar='|', quoting=cs.QUOTE_MINIMAL)
        filewriter.writerow(fieldName)
        filewriter.writerows(finalReturnData)
    print("done")

conjunctive()
disjunctive()
finalReturn()
writeCSV()
