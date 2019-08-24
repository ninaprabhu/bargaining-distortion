#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 22 21:51:43 2017

@author: Nina
"""
import itertools
import numpy as np 
import math 

def all_same(items):
    return all(x == items[0] for x in items) #if all points are the same

def dist(p0, p1):
    return (p0[0] - p1[0])**2 + (p0[1] - p1[1])**2 
    #return abs(((p0[0]-p1[0])**3 + (p0[1]-p1[1])**3)** (1. / 3))

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

def KS(x,y,a):
    # decision reached by agents x and y given alternative a
    # Kalai-Smorodinsky bargaining solution
    def normaldist(p0,p1):
        return math.sqrt((p0[0] - p1[0])**2 + (p0[1] - p1[1])**2)
    if str(x)==str(y) and str(x)==str(a): #if all three points are the same
        return a
    elif str(x)==str(y) or str(x)==str(a): #if two points are the same
        return x
    elif str(a)==str(y): #also if two points are the same
        return y
    elif arecolinear((x,y,a)):
        if dist(x,a)+dist(y,a)==dist(x,y) and dist(x,a)+dist(y,a)==dist(x,y):
            return a
        elif dist(x,y)+dist(a,x)==dist(a,y) and dist(x,y)+dist(a,x)==dist(a,y):
            return x
        else:
            return y
    else: 
        X=np.array(x)
        Y=np.array(y)
        return list(X+(Y-X)*(normaldist(x,a)/(normaldist(y,a)+normaldist(x,a))))

def Omedian(x,y,a):
    if dist(x,a)>dist(y,a):
        o=y
    elif dist(x,a)<dist(y,a):
        o=x
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
    
        
        #if all the elements in the list are the same 
    N=[0]*len(B)
    if all_same(B)==True:
        return 0
        #for some reason the program crashes if all the points are on a line
    elif all_same(collinear(B)):
        z=np.median(B,axis=0)
    else: 
        z=np.mean(B,axis=0)
    #now that we have a point that optimizes social cost, find social cost.
    for u in range(len(B)):
        N[u]= dist(B[u],z) 
        v = sum(N)

    return v

def distortion(C):
    B = list(itertools.permutations(C)) #all permutations of points
    f=[0]*len(B) #empty list to hold distortions
    for n in range(len(B)):
        x = B[n][0] #agent 1
        y = B[n][1] #agent 2
        a = B[n][2] #original alternative
        z = B[n][3] #random agent to take distance from 
        f[n]=dist(KS(x,y,a),z) #two agents bargain, find distance from other agent
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

def Odistortion(C):
    B = list(itertools.permutations(C)) #all permutations of points
    f=[0]*len(B) #empty list to hold distortions
    for n in range(len(B)):
        x = B[n][0] #agent 1
        y = B[n][1] #agent 2
        a = B[n][2] #original alternative
        z = B[n][3] #random agent to take distance from 
        if dist(x,a)==dist(y,a):
            f[n]=0.5*(dist(x,z)+dist(y,z))
        else: 
            f[n]=dist(Omedian(x,y,a),z) #two agents bargain, find distance from other agent
    h = np.mean(f) #mean of all possible arrangements- average SC
    g = optimalsc(C)
    G = g / (len(C)) #average optimal SC
    if g==0: 
        return 1 #if all points are the same
    w = h / G #SC/optimal SC
    return w
