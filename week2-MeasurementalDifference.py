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
#print(*finalData[0])
#print(*ProteinData['HpCD00780223'])
#calculate s2 for markers and average s2
S2MetricData.pop(0)
avg = 0.0
count = 0.0
for Protein in week1Data:
    s2 = float(Protein[1])*float(Protein[2])
    S2MetricData.append([Protein[0],s2])
    avg += s2
    count += 1
averageS2 = avg / count

S2MetricProteinKey = {"test": float}
S2MetricProteinKey.pop("test")
#create subset of markers that have s2 above average
for Protein in S2MetricData:
    if (Protein[1] > averageS2):
        subsetData.append(Protein[0])
        S2MetricProteinKey[Protein[0]] = Protein[1]
#print(S2MetricProteinKey['HpCD00780232'])
subsetData.pop(0)
#print(*subsetData) 


D_r_ = {"test": []}
D_r_.pop("test")
for Protein in subsetData:
    data = ProteinData[Protein]
    tempData = []
    for point in data:
        if (float(point[0]) > 2 ):
            tempData.append(point)
    D_r_[Protein] = tempData
#print(*D_r_['HpCD00780232'])

def Similiarity(Protein1, Protein2):
    dataSet1 = D_r_[Protein1]
    dataSet2 = D_r_[Protein2]
    sharedSamples = 0
    minimumSamples = 1
    for data in dataSet1:
        if data in dataSet2:
            sharedSamples += 1
    if sharedSamples >= minimumSamples:
        return True
    else:
        return False
    #print(dataSet1)
#print (R_r_("HpCD00780232","HpCD00780322"))

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
    #print(p_Y0_D0)
    #print(p_Y1_D1)
    #print()
    return (p_Y1_D1*p_Y0_D0)

def conjunctive():
    R_r_ = {"test": []}
    notR_r_ = {"test": []}
    R_r_.pop("test")
    notR_r_.pop("test")
    for Protein in subsetData:
        tempProteinCorrelatedData = []
        notTempProteinCorrelatedData = []
        for SecondProtein in subsetData:
            if (Similiarity(Protein, SecondProtein) and Protein != SecondProtein):
                dataset1 = ProteinData[Protein]
                dataset2 = ProteinData[SecondProtein]
                #print(dataset1[:50])
                #print(dataset2[:50])
                #print(Protein, "  ", SecondProtein)
                tempSampleData = []
                notTempSampleData = []
                for index in range(len(dataset1)):
                    test = float(dataset1[index][0])
                    #print(*dataset1[index])
                    #print(*dataset2[index])
                    if (float(dataset1[index][0]) > 2 and float(dataset2[index][0]) > 2):
                        tempSampleData.append((1,dataset1[index][2]))
                    else:
                        tempSampleData.append((0, dataset1[index][2]))
                    if (float(dataset1[index][0]) <= 2 and float(dataset2[index][0]) > 2):
                        #print(dataset1[index][0], "   ", dataset2[index][0])
                        notTempSampleData.append((1,dataset1[index][2]))
                    else:
                        notTempSampleData.append((0,dataset1[index][2]))
                #print(tempSampleData[:50])
                #print()
                tempProteinCorrelatedData.append((SecondProtein,tempSampleData))
                notTempProteinCorrelatedData.append((SecondProtein,notTempSampleData))
        if (len(tempProteinCorrelatedData) > 0):
            R_r_[Protein] = tempProteinCorrelatedData
        if (len(notTempProteinCorrelatedData) > 0):
            notR_r_[Protein] = notTempProteinCorrelatedData
    #print(*R_r_)
    for Protein in R_r_:
        eavg = 0
        for index in range(len(R_r_[Protein])):
            dataP_r_r = R_r_[Protein][index]
            datanotP_r_r = notR_r_[Protein][index]
            Pr_r = compute(dataP_r_r[1])
            Pnotr_r = compute(datanotP_r_r[1])
            #print(Pr_r, "   ", Pnotr_r)
            eavg += Pr_r - Pnotr_r
        eavg = eavg/len(R_r_[Protein])
        #print(eavg)
        conjunctiveData[Protein] = eavg

def disjunctive():
    R_r_ = {"test": []}
    notR_r_ = {"test": []}
    R_r_.pop("test")
    notR_r_.pop("test")
    for Protein in subsetData:
        tempProteinCorrelatedData = []
        notTempProteinCorrelatedData = []
        for SecondProtein in subsetData:
            if (Similiarity(Protein, SecondProtein) and Protein != SecondProtein):
                dataset1 = ProteinData[Protein]
                dataset2 = ProteinData[SecondProtein]
                #print(dataset1[:50])
                #print(dataset2[:50])
                print(Protein, "  ", SecondProtein)
                tempSampleData = []
                notTempSampleData = []
                for index in range(len(dataset1)):
                    test = float(dataset1[index][0])
                    #print(*dataset1[index])
                    #print(*dataset2[index])
                    if (float(dataset1[index][0]) > 2 or float(dataset2[index][0]) > 2):
                        tempSampleData.append((1,dataset1[index][2]))
                    else:
                        tempSampleData.append((0, dataset1[index][2]))
                    if (float(dataset1[index][0]) <= 2 or float(dataset2[index][0]) > 2):
                        #print(dataset1[index][0], "   ", dataset2[index][0])
                        notTempSampleData.append((1,dataset1[index][2]))
                    else:
                        notTempSampleData.append((0,dataset1[index][2]))
                #print(tempSampleData[:50])
                #print()
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
            #print(Pr_r, "   ", Pnotr_r)
            eavg += Pr_r - Pnotr_r
        eavg = eavg/len(R_r_[Protein])
        #print(eavg)
        disjunctiveData[Protein] = eavg
    return disjunctiveData

def finalReturn():
    S2Data = S2MetricProteinKey

    # "n/a" means that the specific Protein didnt have a value for conjunctive or disjunctive because it didn't share any common blood markers
    for Protein in S2Data:
        #print(Protein)
        #print(S2MetricProteinKey[Protein])
        #print(*conjunctiveData)
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
#print(*finalReturnData)