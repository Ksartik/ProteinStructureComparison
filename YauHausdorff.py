import sys
from random import uniform
from random import gauss
import math
from numpy.random import permutation as perm
from numpy.random import shuffle 
import pandas as pd
from urllib.request import urlretrieve

protein1 = sys.argv[1]
protein2 = sys.argv[2]

# Converting data from database into dataframe
try :
    f1 = open("PDB/" + protein1 + ".pdb.txt")
except FileNotFoundError :
    urlretrieve("https://files.rcsb.org/view/" + protein1 + ".pdb", "PDB/" + protein1 + ".pdb.txt")
    f1 = open("PDB/" + protein1 + ".pdb.txt")    

# df1 = {"S.No." : [], "Atom" : [], "Atom Type" : [], "AA" : [], "AA no." : [], "Branch" : [], "X" : [], "Y" : [], "Z" : [], "Occupancy" : [], "Temp Factor" : []}
df1 = {"Atom" : [], "Atom Type" : [], "AA" : [], "AA no." : [], "Branch" : [], "X" : [], "Y" : [], "Z" : [], "Occupancy" : [], "Temp Factor" : []}
for line in f1.readlines() :
    # data = list(filter(lambda x : x != '', line.split(" ")))
    # if (data[0] == 'ATOM') :
        # df1["S.No."].append(data[1])
        # df1["Atom Type"].append(data[2])
        # df1["AA"].append(data[3])
        # df1["Branch"].append(data[4])
        # df1["AA no."].append(int(data[5]))
        # df1["X"].append(float(data[6]))
        # df1["Y"].append(float(data[7]))
        # df1["Z"].append(float(data[8]))
        # df1["Occupancy"].append(float(data[9]))
        # df1["Temp Factor"].append(float(data[10]))
        # df1["Atom"].append(data[11])
    if (line[0:4] == "ATOM") :
        df1["Atom Type"].append(line[13:17].strip())
        df1["AA"].append(line[17:21].strip())
        df1["Branch"].append(line[21:22])
        df1["AA no."].append(int(line[22:26].strip()))
        df1["X"].append(float(line[26:38].strip()))
        df1["Y"].append(float(line[38:46].strip()))
        df1["Z"].append(float(line[46:54].strip()))
        df1["Occupancy"].append(float(line[54:60].strip()))
        df1["Temp Factor"].append(float(line[60:66].strip()))
        df1["Atom"].append(line[77:].strip())

df1 = pd.DataFrame(df1)


#same with structure 2
try :
    f2 = open("PDB/" + protein2 + ".pdb.txt")
except FileNotFoundError :
    urlretrieve("https://files.rcsb.org/view/" + protein2 + ".pdb", "PDB/" + protein2 + ".pdb.txt")
    f2 = open("PDB/" + protein2 + ".pdb.txt")    

# df2 = {"S.No." : [], "Atom" : [], "Atom Type" : [], "AA" : [], "AA no." : [], "Branch" : [], "X" : [], "Y" : [], "Z" : [], "Occupancy" : [], "Temp Factor" : []}
df2 = {"Atom" : [], "Atom Type" : [], "AA" : [], "AA no." : [], "Branch" : [], "X" : [], "Y" : [], "Z" : [], "Occupancy" : [], "Temp Factor" : []}
for line in f2.readlines() :
    # data = list(filter(lambda x : x != '', line.split(" ")))
    # if (data[0] == 'ATOM') :
        # df2["S.No."].append(data[1])
        # df2["Atom Type"].append(data[2])
        # df2["AA"].append(data[3])
        # df2["Branch"].append(data[4])
        # df2["AA no."].append(int(data[5]))
        # df2["X"].append(float(data[6]))
        # df2["Y"].append(float(data[7]))
        # df2["Z"].append(float(data[8]))
        # df2["Occupancy"].append(float(data[9]))
        # df2["Temp Factor"].append(float(data[10]))
        # df2["Atom"].append(data[11])
    if (line[0:4] == "ATOM") :
        df2["Atom Type"].append(line[13:17].strip())
        df2["AA"].append(line[17:21].strip())
        df2["Branch"].append(line[21:22])
        df2["AA no."].append(int(line[22:26].strip()))
        df2["X"].append(float(line[26:38].strip()))
        df2["Y"].append(float(line[38:46].strip()))
        df2["Z"].append(float(line[46:54].strip()))
        df2["Occupancy"].append(float(line[54:60].strip()))
        df2["Temp Factor"].append(float(line[60:66].strip()))
        df2["Atom"].append(line[77:].strip())

df2 = pd.DataFrame(df2)


## generating points from protein data frames

structure1 = [] 
structure2 = []
for i in range(len(df1.index)) :
    structure1.append((df1.iloc[i]["X"], df1.iloc[i]["Y"], df1.iloc[i]["Z"]))
    
for i in range(len(df2.index)) :
    structure2.append((df2.iloc[i]["X"], df2.iloc[i]["Y"], df2.iloc[i]["Z"]))

## checking for random points --
# random point generation

def uniformArray (lb, ub, num) :
    A = []
    for i in range(num) :
        A.append(uniform(lb, ub))
    return(A)

def uniformPoints (lb, ub, num) :
    Ax = uniformArray(lb, ub, num)
    Ay = uniformArray(lb, ub, num)
    Az = uniformArray(lb, ub, num)
    return (list(map(lambda x,y,z : (x,y,z), Ax, Ay, Az)))

A = uniformPoints(-1000, 1000, 1000)
B = uniformPoints(-1000, 1000, 100)

## end of random points generation 

# algorithm for Yau Hausdorff distance calculation starts - 

def Bsearch (x, A, acc) :
    large = len(A) - 1
    small = 0
    while (small <= large) :
        mid = (small + large) // 2
        amid = math.floor(A[mid] * (10**acc))
        accx = math.floor(x * (10**acc))
        if (amid < accx) :
            small = mid + 1
        elif (amid > accx) :
            large = mid - 1
        else :
            return (mid, mid)
    return (small - 1, small)

# def f (A, B, t) :
#     if (-A[0] >= B[-1]) :
#         if (t <= (B[-1] - A[0])/2) :
#             return (- t - A[0])
#         else :
#             return (t - B[-1])
#     else :
#         if (t <= (B[-1] - A[0])/2) :
#             return (B[-1] - t)
#         else :
#             return (A[0] + t)

def minHausdorff (A, B, acc = 5) :
    A.sort()
    B.sort()
    if (A[-1] - A[0] < B[-1] - B[0]) :
        temp = B
        B = A
        A = temp
    A = list (map(lambda x : x - A[-1], A))
    B = list (map(lambda x : x - B[0], B))
    s1 = []
    s2 = []
    for i in range(len(A)) :
        sresult1 = Bsearch(B[-1] + A[i], B, acc)
        sresult2 = Bsearch(-A[0] + A[i], B, acc)
        if ((sresult1 == sresult2) and (sresult1[0]> -1) and (sresult1[1] < len(B))):
            s1.append((B[sresult1[0]] - A[i], B[sresult1[1]] - A[i]))
    Aminus = list (map(lambda x : -x, A))
    Aminus = Aminus[::-1]
    for j in range(len(B)) :
        sresult1 = Bsearch(B[-1] - B[j], Aminus, acc)
        sresult2 = Bsearch(-A[0] - B[j], Aminus, acc)
        if ((sresult1 == sresult2) and (sresult1[0] > -1) and (sresult1[1] < len(Aminus))) :
            s2.append((Aminus[sresult1[0]] + B[j], Aminus[sresult1[1]] + B[j]))
    s = list(set(s1) | set(s2))
    s.sort(key = lambda x : x[0])
    for i in range(len(s)) :
        j = i + 1
        while ((j < len(s)) and (s[j][0] < s[i][1])):
            if (s[j][1] < s[i][1]) :
                s.pop(j)
            j = j + 1
    mindist = abs(-(B[-1] + A[0])/2)
    for inter in s :
        if (B[-1] <= -A[0]) :
            d = min ([(inter[0] - A[0])/2, (inter[1] + B[-1])/2])
        else :
            d = min ([(inter[0] + B[-1])/2, (inter[1] - A[0])/2])
        if (d < mindist) :
            mindist = d
    return mindist

def pxAtheta (A, theta) :
    return (list (map(lambda p : (p[0]*math.cos(theta) - p[1]*math.sin(theta)), A)))

def pxBphi (B, phi) :
    return (list (map(lambda p : (p[0]*math.cos(phi) + p[2]*math.sin(phi)), B)))

def ebHausdorff (A, B, metric) :
    cmax = 0
    Ar = list(perm(A))
    Br = list(perm(B))
    for x in Ar :
        cmin = math.inf
        for y in Br :
            d = metric(x, y) 
            if (d < cmax) :
                cmin = 0
                break
            cmin = min([cmin, d])
        cmax = max([cmax, cmin])
    return cmax

def yauHausdorff(A, B) :
    theta_set = [pxAtheta(A, 0.0)]
    phi_set = [pxBphi(B, 0.0)]
    N = 50
    for i in range(N) :
        theta_set.append(pxAtheta(A, uniform(2*math.pi*i/N, 2*math.pi*(i+1)/N)))
        phi_set.append(pxBphi(B, uniform(2*math.pi*i/N, 2*math.pi*(i+1)/N)))
        # theta_set.append(pxAtheta(A, math.pi*(2*i+1)/N))
        # phi_set.append(pxBphi(B, math.pi*(2*i+1)/N))
    return (ebHausdorff(theta_set, phi_set, minHausdorff))
    # return (minHausdorff(theta_set[1], phi_set[1]))


# random point checking 
print(yauHausdorff(structure1, structure1))
# print("Yau Hausdorff for A vs A : ", yauHausdorff(A,A))
# Aold = A
# shuffle(A)
# print("Yau Hausdorff for A vs shuffled A : ", yauHausdorff(A,A))
# print("Yau Hausdorff for A vs B : ", yauHausdorff(Aold, B))
# print("Yau Hausdorff for shuffled A vs B : ", yauHausdorff(A, B))
