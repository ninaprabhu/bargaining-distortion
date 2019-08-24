 h#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  7 12:39:02 2017

@author: Nina
"""
import numpy as np
import itertools

def all_same(items):
    return all(x == items[0] for x in items) #if all points are the same

def dist(p0, p1):
    import math
    return math.sqrt((p0[0] - p1[0])**2 + (p0[1] - p1[1])**2) 

def slope(p0,p1):
    if p0[0]-p1[0]==0:
        return "vertical"
    else:
        return (p0[1]-p1[1])/(p0[0]-p1[0]) 

def mymedian(x,y,a):
    # decision reached by agents x and y given alternative a
    # Nash bargaining solution
    if str(x)==str(y) and str(x)==str(a): #if all three points are the same
        o = a
    elif str(x)==str(y) or str(x)==str(a): #if two points are the same
        o=x
    elif str(a)==str(y): #also if two points are the same
        o=y
    else:
        X = np.array(x)
        Y = np.array(y)
        k = list((Y-X)/dist(x, y)) #unit vector between x and y
        t = (dist(x,y)+dist(a,x)-dist(a,y))/2 #ratio to multiply unit vector by to get Nash equilibrium
        o = [x[0]+k[0]*t,x[1]+k[1]*t] #finds right point
    return o

def arecolinear(points):
    #finding out if three points in the set are on a line 
    if str(points[0])==str(points[1]) or str(points[0])==str(points[2]) or str(points[2])==str(points[1]):
        return True
    else: 
        slope1 = slope(points[0],points[1])
        slope2 = slope(points[1],points[2])

        if slope1 == slope2:
            return True
        else:
            return False

def collinear(pointlist):
    #taking 
    D=list(itertools.combinations(pointlist,3))
    E=[0]*len(D)
    for n in range(len(D)):
        if arecolinear(D[n]):
            E[n]='Collinear'
        else: 
            E[n]='Not collinear'
    return E

def optimalsc(B):
    # The social cost of the optimal decision
    
        
        #if all the elements in the matrix are the same 
    N=[0]*len(B)
    if all_same(B)==True:
        return 0
        #for some reason the program crashes if all the points are on a line
    elif all_same(collinear(B)):
        z=np.median(B,axis=0)
    else: 
        #iterating with Weiszfeld's algorithm
        s=[[0,0]]*len(B)
        t=[0]*len(B)
        for u in range(len(B)):
            s[u] = [(B[u][0] / (dist(B[u], (.3,.3)))),(B[u][1] / dist(B[u], (.3,.3)))]
            S    = np.sum(s,axis=0)
            t[u] = 1 / dist(B[u],(.3 ,.3))
            T    = sum(t)
        a   =[0]*2
        a[0]= S[0]/T
        a[1]= S[1]/T
        if tuple(a) in B:
            z=a
        else: 
            s=[(0,0)]*len(B)
            t=[0]*len(B)
            for u in range(len(B)):
                s[u] = [B[u][0] / (dist(B[u], a)),(B[u][1] / (dist(B[u], a)))]
                S    = np.sum(s,axis=0)
                t[u] = 1 / dist(B[u],a)
                T    = sum(t)
            b   =[0]*2
            b[0]= S[0]/T
            b[1]= S[1]/T
            if tuple(b) in B:
                z=b
            else: 
                s=[(0,0)]*len(B)
                t=[0]*len(B)
                for u in range(len(B)):
                    s[u] = (B[u][0] / (dist(B[u], b)),(B[u][1] / (dist(B[u], b))))
                    S    = np.sum(s,axis=0)
                    t[u] = 1 / dist(B[u],b)
                    T    = sum(t)
                c   =[0]*2
                c[0]= S[0]/T
                c[1]= S[1]/T
                if tuple(c) in B:
                    z=c
                else: 
                    s=[(0,0)]*len(B)
                    t=[0]*len(B)
                    for u in range(len(B)):
                        s[u] = (B[u][0] / (dist(B[u], c)),(B[u][1] / (dist(B[u], c))))
                        S    = np.sum(s,axis=0)
                        t[u] = 1 / dist(B[u],c)
                        T    = sum(t)
                    d   =[0]*2
                    d[0]= S[0]/T
                    d[1]= S[1]/T
                    if tuple(d) in B:
                        z=d
                    else: 
                        s=[(0,0)]*len(B)
                        t=[0]*len(B)
                        for u in range(len(B)):
                            s[u] = (B[u][0] / (dist(B[u], d)),(B[u][1] / (dist(B[u], d))))
                            S    = np.sum(s,axis=0)
                            t[u] = 1 / dist(B[u],d)
                            T    = sum(t)
                        e   =[0]*2
                        e[0]= S[0]/T
                        e[1]= S[1]/T
                        if tuple(e) in B:
                            z=e
                        else: 
                            s=[(0,0)]*len(B)
                            t=[0]*len(B)
                            for u in range(len(B)):
                                s[u] = (B[u][0] / (dist(B[u], e)),(B[u][1] / (dist(B[u], e))))
                                S    = np.sum(s,axis=0)
                                t[u] = 1 / dist(B[u],e)
                                T    = sum(t)
                            f   =[0]*2
                            f[0]= S[0]/T
                            f[1]= S[1]/T
                            if tuple(f) in B:
                                z=f
                            else: 
                                s=[(0,0)]*len(B)
                                t=[0]*len(B)
                                for u in range(len(B)):
                                    s[u] = (B[u][0] / (dist(B[u], f)),(B[u][1] / (dist(B[u], f))))
                                    S    = np.sum(s,axis=0)
                                    t[u] = 1 / dist(B[u],f)
                                    T    = sum(t)
                                g   =[0]*2
                                g[0]= S[0]/T
                                g[1]= S[1]/T
                                if tuple(g) in B: 
                                    z=g
                                else:
                                    s=[(0,0)]*len(B)
                                    t=[0]*len(B)
                                    for u in range(len(B)):
                                        s[u] = (B[u][0] / (dist(B[u], g)),(B[u][1] / (dist(B[u], g))))
                                        S    = np.sum(s,axis=0)
                                        t[u] = 1 / dist(B[u],g)
                                        T    = sum(t)
                                    h   =[0]*2
                                    h[0]= S[0]/T
                                    h[1]= S[1]/T
                                    if tuple(h) in B:
                                        z=h
                                    else:
                                        s=[(0,0)]*len(B)
                                        t=[0]*len(B)
                                        for u in range(len(B)):
                                            s[u] = (B[u][0] / (dist(B[u], h)),(B[u][1] / (dist(B[u], h))))
                                            S    = np.sum(s,axis=0)
                                            t[u] = 1 / dist(B[u],h)
                                            T    = sum(t)
                                        j  =[0]*2
                                        j[0]= S[0]/T
                                        j[1]= S[1]/T
                                        if tuple(j) in B:
                                            z=j
                                        else: 
                                            s=[(0,0)]*len(B)
                                            t=[0]*len(B)
                                            for u in range(len(B)):
                                                s[u] = (B[u][0] / (dist(B[u], j)),(B[u][1] / (dist(B[u], j))))
                                                S    = np.sum(s,axis=0)
                                                t[u] = 1 / dist(B[u],j)
                                                T    = sum(t)
                                            k   =[0]*2
                                            k[0]= S[0]/T
                                            k[1]= S[1]/T
                                            if tuple(k) in B:
                                                z=k
                                            else:
                                                s=[(0,0)]*len(B)
                                                t=[0]*len(B)
                                                for u in range(len(B)):
                                                    s[u] = (B[u][0] / (dist(B[u], k)),(B[u][1] / (dist(B[u], k))))
                                                    S    = np.sum(s,axis=0)
                                                    t[u] = 1 / dist(B[u],k)
                                                    T    = sum(t)
                                                l   =[0]*2
                                                l[0]= S[0]/T
                                                l[1]= S[1]/T
                                                if tuple(l) in B:
                                                    z=l
                                                else: 
                                                    s=[(0,0)]*len(B)
                                                    t=[0]*len(B)
                                                    for u in range(len(B)):
                                                        s[u] = (B[u][0] / (dist(B[u], l)),(B[u][1] / (dist(B[u], l))))
                                                        S    = np.sum(s,axis=0)
                                                        t[u] = 1 / dist(B[u],l)
                                                        T    = sum(t)
                                                    m   =[0]*2
                                                    m[0]= S[0]/T
                                                    m[1]= S[1]/T
                                                    if tuple(m) in B:
                                                        z=m
                                                    else: 
                                                        s=[(0,0)]*len(B)
                                                        t=[0]*len(B)
                                                        for u in range(len(B)):
                                                            s[u] = (B[u][0] / (dist(B[u], m)),(B[u][1] / (dist(B[u], m))))
                                                            S    = np.sum(s,axis=0)
                                                            t[u] = 1 / dist(B[u],m)
                                                            T    = sum(t)
                                                        n   =[0]*2
                                                        n[0]= S[0]/T
                                                        n[1]= S[1]/T
                                                        if tuple(n) in B:
                                                            z=n
                                                        else: 
                                                            s=[(0,0)]*len(B)
                                                            t=[0]*len(B)
                                                            for u in range(len(B)):
                                                                s[u] = (B[u][0] / (dist(B[u], n)),(B[u][1] / (dist(B[u], n))))
                                                                S    = np.sum(s,axis=0)
                                                                t[u] = 1 / dist(B[u],n)
                                                                T    = sum(t)
                                                            o   =[0]*2
                                                            o[0]= S[0]/T
                                                            o[1]= S[1]/T
                                                            if tuple(o) in B:
                                                                z=o
                                                            else: 
                                                                s=[(0,0)]*len(B)
                                                                t=[0]*len(B)
                                                                for u in range(len(B)):
                                                                    s[u] = (B[u][0] / (dist(B[u], o)),(B[u][1] / (dist(B[u], o))))
                                                                    S    = np.sum(s,axis=0)
                                                                    t[u] = 1 / dist(B[u],o)
                                                                    T    = sum(t)
                                                                p   =[0]*2
                                                                p[0]= S[0]/T
                                                                p[1]= S[1]/T
                                                            if tuple(p) in B:
                                                                z=p
                                                            else:
                                                                s=[(0,0)]*len(B)
                                                                t=[0]*len(B)
                                                                for u in range(len(B)):
                                                                    s[u] = (B[u][0] / (dist(B[u], p)),(B[u][1] / (dist(B[u], p))))
                                                                    S    = np.sum(s,axis=0)
                                                                    t[u] = 1 / dist(B[u],p)
                                                                    T    = sum(t)
                                                                q   =[0]*2
                                                                q[0]= S[0]/T
                                                                q[1]= S[1]/T
                                                                if tuple(q) in B:
                                                                    z=q
                                                                else: 
                                                                    s=[(0,0)]*len(B)
                                                                    t=[0]*len(B)
                                                                    for u in range(len(B)):
                                                                        s[u] = (B[u][0] / (dist(B[u], q)),(B[u][1] / (dist(B[u], q))))
                                                                        S    = np.sum(s,axis=0)
                                                                        t[u] = 1 / dist(B[u],q)
                                                                        T    = sum(t)
                                                                    r   =[0]*2
                                                                    r[0]= S[0]/T
                                                                    r[1]= S[1]/T
                                                                    if tuple(r) in B:
                                                                        z=r
                                                                    else: 
                                                                        s=[(0,0)]*len(B)
                                                                        t=[0]*len(B)
                                                                        for u in range(len(B)):
                                                                            s[u] = (B[u][0] / (dist(B[u], r)),(B[u][1] / (dist(B[u], r))))
                                                                            S    = np.sum(s,axis=0)
                                                                            t[u] = 1 / dist(B[u],r)
                                                                            T    = sum(t)
                                                                        x   =[0]*2
                                                                        x[0]= S[0]/T
                                                                        x[1]= S[1]/T
                                                                        if tuple(x) in B:
                                                                            z=x
                                                                        else: 
                                                                            for u in range(len(B)):
                                                                                s[u] = (B[u][0] / (dist(B[u], x)),(B[u][1] / (dist(B[u], x))))
                                                                                S    = np.sum(s,axis=0)
                                                                                t[u] = 1 / dist(B[u],x)
                                                                                T    = sum(t)
                                                                            y   =[0]*2
                                                                            y[0]= S[0]/T
                                                                            y[1]= S[1]/T
                                                                            if tuple(y) in B:
                                                                                z=y
                                                                            else: 
                                                                                s=[(0,0)]*len(B)
                                                                                t=[0]*len(B)
                                                                                for u in range(len(B)):
                                                                                    s[u] = (B[u][0] / (dist(B[u], y)),(B[u][1] / (dist(B[u], y))))
                                                                                    S    = np.sum(s,axis=0)
                                                                                    t[u] = 1 / dist(B[u],y)
                                                                                    T    = sum(t)
                                                                                finalminimizer=[0]*2
                                                                                finalminimizer[0]= S[0]/T
                                                                                finalminimizer[1]= S[1]/T
                                                                                z=finalminimizer
    #now that we have a point that optimizes social cost, find social cost.
    for u in range(len(B)):
        N[u]= dist(B[u],z) 
        v = sum(N)

    return v

def singledistortion(points): #from when I tried to do something that applies function to every row in array
    return dist(mymedian(points[0], points[1], points[2]), points[3])
    

def distortion(C):
    B = list(itertools.permutations(C)) #all permutations of points
    f=[0]*len(B) #empty list to hold distortions
    for n in range(len(B)):
        x = B[n][0] #agent 1
        y = B[n][1] #agent 2
        a = B[n][2] #original alternative
        z = B[n][3] #random agent to take distance from 
        f[n]=dist(mymedian(x,y,a),z) #two agents bargain, find distance from other agent
    h = np.mean(f) #mean of all possible arrangements- average SC
    g = optimalsc(C)
    G = g / (len(C)) #average optimal SC
    if g==0: 
        return 1 #if all points are the same
    w = h / G #SC/optimal SC
    return w


def grid(a,b):
    #b points are chosen from an a x a grid
    C = range(a+1) #numbers to choose from
    B = list(itertools.product(C,repeat=2)) #all possible points
    middlepoints = list(itertools.combinations_with_replacement(B,b-3)) #chooses b-3 points
    left = [None]*(a+1) #empty 
    top  = [None]*(a+1) #empty
    bottom = [None]*(a+1) #empty
    for n in range(a+1):
        left[n]  = (0,n) #choose a point on y axis
        bottom[n]= (n,0) #bottom of grid- x axis
        top[n]   = (n,a) #top of grid
        
    combo = (left,top,bottom) 
    sidepoints   = list(itertools.product(*combo)) #all combinations of side points
    middleandside= [middlepoints,sidepoints] 
    allpoints=list(itertools.product(*middleandside)) #all possible sets of b points
    allpoints2=[None]*len(allpoints)      
    for n in range(len(allpoints)):
        allpoints2[n]=allpoints[n][0]+allpoints[n][1] #allpoints didn't combine them, so this one does
    f=[0]*2
    for n in range(len(allpoints2)):
        s=distortion(allpoints2[n]) #assigns distortion to a variable
        if s>f[1]: 
            f=[allpoints2[n],s] #replaces current distortion with new if bigger
    return f #return biggest distortion and points that make it
