import sys
from random import uniform
import math
from numpy.random import permutation as perm
from numpy.random import shuffle 
# protein1 = sys.argv[2]
# protein2 = sys.argv[3]
# # structure1 = {"atoms" : [], "coordinates" : []}
# structure1 = {}
# with open ("Data/" + protein1 + ".txt") as file:
#     for line in file.readlines() :
#         data = line.split("    ")
        # structure1["atoms"].append(data[0]) 
        # structure1["coordinates"].append([data[1], data[2], data[3]]) #parse integer
        

#same with structure 2



## checking for random points --

# random point generation

def uniformArray (lb, ub, length) :
    A = []
    for i in range(length) :
        A.append(uniform(lb, ub))
    return(A)

def uniformPoints (lb, ub, length) :
    Ax = uniformArray(lb, ub, length)
    Ay = uniformArray(lb, ub, length)
    Az = uniformArray(lb, ub, length)
    return (list(map(lambda x,y,z : (x,y,z), Ax, Ay, Az)))

A = uniformPoints(-1000, 1000, 1000)
B = uniformPoints(-1000, 1000, 1000)

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
    return (list (map(lambda p : (p[0]*math.cos(theta) + p[1]*math.sin(theta)), A)))

def pxBphi (B, phi) :
    return (list (map(lambda p : (p[0]*math.cos(phi) - p[2]*math.sin(phi)), B)))

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
        phi_set.append(pxAtheta(B, uniform(2*math.pi*i/N, 2*math.pi*(i+1)/N)))
    return (ebHausdorff(theta_set, phi_set, minHausdorff))


# random point checking 

print("Yau Hausdorff for A vs A : ", yauHausdorff(A,A))
Aold = A
shuffle(A)
print("Yau Hausdorff for A vs shuffled A : ", yauHausdorff(A,A))
print("Yau Hausdorff for A vs B : ", yauHausdorff(Aold, B))
print("Yau Hausdorff for shuffled A vs B : ", yauHausdorff(A, B))
