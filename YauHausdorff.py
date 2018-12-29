import sys
from random import uniform
from random import gauss
import math
import numpy as np
from numpy.random import permutation as perm
from numpy.random import shuffle 
from numpy.random import seed
import pandas as pd
from urllib.request import urlretrieve
import randrot
import matplotlib.pyplot as plt

protein1 = sys.argv[1]
protein2 = sys.argv[2]

# Converting data from database into dataframe

def readPDB (protein) : 
    try :
        f1 = open("PDB/" + protein + ".pdb.txt")
    except FileNotFoundError :
        urlretrieve("https://files.rcsb.org/view/" + protein + ".pdb", "PDB/" + protein + ".pdb.txt")
        f1 = open("PDB/" + protein + ".pdb.txt")
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
    return (pd.DataFrame(df1))

df1 = readPDB(protein1)
df2 = readPDB(protein2)

## generating points from protein data frames

structure1 = [] 
structure2 = []
for i in range(len(df1.index)) :
    structure1.append((df1.iloc[i]["X"], df1.iloc[i]["Y"], df1.iloc[i]["Z"]))
    
for i in range(len(df2.index)) :
    structure2.append((df2.iloc[i]["X"], df2.iloc[i]["Y"], df2.iloc[i]["Z"]))

## checking for random points --
# random point generation

seed(0)
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

A = uniformPoints(-100, 100, 100)
B = uniformPoints(-100, 100, 100)

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

def generalRot(A, theta3) :
    rx = [[1,0,0], [0, math.cos(theta3[0]), -math.sin(theta3[0])], [0, math.sin(theta3[0]), math.cos(theta3[0])]]
    ry = [[math.cos(theta3[1]), 0, math.sin(theta3[1])], [0, 1, 0], [-math.sin(theta3[1]), 0, math.cos(theta3[1])]]
    rz = [[math.cos(theta3[2]), math.sin(theta3[2]), 0], [-math.sin(theta3[2]), math.cos(theta3[2]), 0], [0, 0, 1]]
    return (np.dot(A, np.dot(rx, np.dot(ry, rz))))

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

def bruteForceYH (A, B) :
    amax = 0
    for pa in A :
        bmin = math.inf
        for pb in B :
            d = minHausdorff(pa, pb)
            bmin = min([bmin, d])
        amax = max([amax, bmin])
    return amax

def centroid(X):
    C = X.mean(axis=0)
    return C

def kabsch(P, Q):
    # Computation of the covariance matrix
    P = (P - centroid(P))
    Q = (Q - centroid(Q))
    C = np.dot(np.transpose(P), Q)

    # Computation of the optimal rotation matrix
    # This can be done using singular value decomposition (SVD)
    # Getting the sign of the det(V)*(W) to decide
    # whether we need to correct our rotation matrix to ensure a
    # right-handed coordinate system.
    # And finally calculating the optimal rotation matrix U
    V, S, W = np.linalg.svd(C)
    d = (np.linalg.det(V) * np.linalg.det(W)) < 0.0

    if d:
        S[-1] = -S[-1]
        V[:, -1] = -V[:, -1]

    # Create Rotation matrix U
    U = np.dot(V, W)

    return U

def rotate(P, R) :
    return np.dot(P, R)

def bruteForceHDD3 (A, B) :
    dist_mat = np.zeros(shape = (len(A), len(B)))
    # Directed HDD from A to B
    amax = 0
    for i in range(0, len(A)) :
        bmin = math.inf
        for j in range(0, len(B)) :
            d = dist_mat[i][j] = minHausdorff(A[i], B[j])
            if (d < bmin) :
                bmin = d
        if (bmin > amax) :
            amax = bmin
    # Directed HDD from B to A
    bmax = 0
    for i in range(0, len(B)) :
        amin = math.inf
        for j in range(0, len(A)) :
            d = dist_mat[j][i]
            if (d < amin) :
                amin = d
        if (amin > bmax) :
            bmax = amin
    return (max ([amax, bmax]))

def ebHDD3 (A, B) :
    dist_mat = np.full((len(A), len(B)), fill_value = np.nan)

    # directed hdd from A to B
    cmax = 0
    Ar = list(perm(A))
    Br = list(perm(B))
    for i in range(len(Ar)) :
        cmin = math.inf
        for j in range(len(Br)) :
            d = dist_mat[i][j] = minHausdorff(Ar[i], Br[j])
            if (d < cmax) :
                cmin = 0
                break
            cmin = min([cmin, d])
        cmax = max([cmax, cmin])
    
    # directed hdd from B to A
    dmax = 0
    for i in range(len(Ar)) :
        dmin = math.inf
        for j in range(len(Br)) :
            if (np.isnan(dist_mat[i][j])) :
                d = dist_mat[i][j] = minHausdorff(Ar[i], Br[j])
            else :
                d = dist_mat[i][j]
            if (d < dmax) :
                dmin = 0
                break
            dmin = min([dmin, d])
        dmax = max([dmax, dmin])
    
    return (max([dmax, cmax]))

def yauHausdorff(A, B, s = 0) :
    # theta_set = [pxAtheta(A, 0.0)]
    # phi_set = [pxBphi(B, 0.0)]
    theta_set = []
    phi_set = []
    N = 50
    # for i in range(N) :
    #     theta1 = uniform (2*math.pi*i/N, 2*math.pi*(i+1)/N)
    #     Aarr = np.array(A)
    #     Barr = np.array(B)
    #     # R0 = kabsch (Aarr, Barr)
    #     A1 = rotate(A, np.array([[math.cos(theta1), math.sin(theta1), 0], [-math.sin(theta1), math.cos(theta1), 0], [0,0,1]]))
    #     # A1 = rotate(Aarr, R0)
    #     # R is the rotation matrix that minimizes RMSD b/w A1 and B 
    #     R = kabsch(Barr, A1)
    #     theta_set.append(A1[:,0].tolist())
    #     phi_set.append(rotate(Barr, R)[:,0].tolist())
    #     # theta_set.append(pxAtheta(A, theta1))
    #     # phi_set.append(pxBphi(B, uniform(2*math.pi*i/N, 2*math.pi*(i+1)/N)))
    #     # x = uniform(2*math.pi*i/N, 2*math.pi*(i+1)/N)
    #     # theta_set.append(pxAtheta(A, x))
    #     # phi_set.append(pxBphi(B, x))
    # return (ebHausdorff(theta_set, phi_set, minHausdorff))
    # return (minHausdorff(theta_set[1], phi_set[1]))
    seed(s)
    for i in range(N) :
        theta_set.append(np.asarray(randrot.generate_3d()))
        phi_set.append(np.asarray(randrot.generate_3d()))
    Xa = list(map(lambda theta : rotate(A, theta)[:, 0].tolist(), theta_set))
    Xb = list(map(lambda phi : rotate(B, phi)[:, 0].tolist(), phi_set))
    # return (max([ebHausdorff(Xa, Xb, minHausdorff), ebHausdorff(Xb, Xa, minHausdorff)]))
    # return (max([bruteForceYH(Xa, Xb), bruteForceYH(Xb, Xa)]))
    return bruteForceHDD3(Xa, Xb)

def minCell (A, B, N, itheta = 0, iphi = 0, ftheta = 2*math.pi, fphi = 2*math.pi) :
    theta_set = []
    phi_set = []
    for i in range(N) :
        theta_set.append(itheta + (ftheta - itheta)*i/N)
        phi_set.append(iphi + (fphi - iphi)*i/N)
    # Xa = list(map(lambda theta : rotate(A, theta)[:, 0].tolist(), theta_set))
    # Xb = list(map(lambda phi : rotate(B, phi)[:, 0].tolist(), phi_set))
    Xa = list(map(lambda theta : pxAtheta(A, theta), theta_set))
    Xb = list(map(lambda phi : pxBphi(B, phi), phi_set))
    mindis, mintheta, minphi = (math.inf, 0, 0)
    for i in range(len(Xa)) :
        for j in range(len(Xb)) :
            dis = minHausdorff(Xa[i], Xb[j])
            if (dis < mindis) :
                mindis = dis
                mintheta = theta_set[i]
                minphi = phi_set[j]
    return (mintheta, minphi)

def yhMittal (A, B) :
    # N = 18
    N = 9
    itheta = 0
    iphi = 0
    ftheta = 2*math.pi
    fphi = 2*math.pi
    mtheta, mphi = (math.pi, math.pi)
    while ((abs(ftheta - itheta) < math.pi/180) and (abs(fphi - iphi) < math.pi/180)) :
        mtheta, mphi = minCell(A, B, N, itheta, iphi, ftheta, fphi)
        tdiff = (ftheta - itheta)/N
        itheta = mtheta - tdiff
        ftheta = mtheta + tdiff
        pdiff = (fphi - iphi)/N
        iphi = mphi - pdiff
        fphi = mphi + pdiff
    k = 50
    theta_set = list(map(lambda x : mtheta + x*(ftheta - itheta)/k, range(-math.floor(k/2), 1 + math.floor(k/2))))
    phi_set = list(map(lambda x : mphi + x*(fphi - iphi)/k, range(-math.floor(k/2), 1 + math.floor(k/2))))
    Xa = list(map(lambda theta : pxAtheta(A, theta), theta_set))
    Xb = list(map(lambda phi : pxBphi(B, phi), phi_set))
    return bruteForceHDD3(Xa, Xb)

        
# print(yauHausdorff(structure1, structure2, 150))
print (yhMittal(structure1, structure2))

'''
            Results for yhMittal (with N = 9 rather than 18 - one can get better and better) :: 

        1R93    1WFZ    2D2A    2JNV
1R93    5.488   7.410   7.942   6.768
1WFZ    7.410   6.513   7.847   18.41
2D2A    7.942   7.847   17.09   15.84
2JNV    6.768   18.41   3.059   4.319

'''



# random point checking 
# X = range(10)
# Y = list(map(lambda x : yauHausdorff(A, A, x), X))
# # plt.plot(X, Y)
# # plt.show()
# ymin = Y[0]
# minx = 0
# for i in range(1, len(X)):
#     y = Y[i]
#     x = X[i]
#     if (y < ymin) :
#         minx = x
#         ymin = y

# print (minx)

# print("Yau Hausdorff for A vs A : ", yauHausdorff(A,A))
# Aold = A
# shuffle(A)
# print("Yau Hausdorff for A vs shuffled A : ", yauHausdorff(A,A))
# print("Yau Hausdorff for A vs B : ", yauHausdorff(Aold, B))
# print("Yau Hausdorff for shuffled A vs B : ", yauHausdorff(A, B))