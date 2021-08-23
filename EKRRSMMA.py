from numpy import *
import csv
#import math
#import pandas as pd
import random 
import time
#olderr = seterr(all='ignore')
nm = 541 
ns = 831 
nc = 664 
r = 0.8 
M = 20 
g= 0.7 
adjustparameter =0.5 # adjustable parameter
regularpara = 1 # regularization parameter

SMS = loadtxt(r'SM similarity matrix.txt') 
ms = loadtxt(r'miRNA similarity matrix.txt')  

ConnectDate = loadtxt(r'known SM-miRNA associations.txt',dtype=int)-1 
SM_namenumber = loadtxt(r'SM number.txt',dtype=bytes).astype(str)
miRNA_namenumber= loadtxt(r'miRNA number.txt',dtype=bytes).astype(str)

def svdprocess(s,row):
   
    U,Sigma,V=linalg.svd(s) 
    k=int(row*0.2)  
    
    SigmaDimred=Sigma.reshape(row,1)[0:k]
    UDimred=U[:,:k]
    SigmaDimred_matrix = zeros((k,k)) 
    for i in range(k):
        SigmaDimred_matrix[i,i] =Sigma.reshape(row,1)[i]
    s_Dimred=dot(UDimred,SigmaDimred_matrix) 
    return s_Dimred


A = zeros((ns,nm),dtype=float)    
for i in range(nc):
    A[ConnectDate[i,0], ConnectDate[i,1]] = 1

start=time.time()


rowsum = sum(A,axis=1)
colsum = sum(A,axis=0) 
row0=argwhere(rowsum==0).ravel() 
col0=argwhere(colsum==0).ravel()   

As_WNN=A.copy()
Am_WNN=A.copy()    
for i in row0:
    for j in range(ns):
        As_WNN[i] = As_WNN[i] + (g**j)*A[argsort(-SMS[i])[j]] #WNN infers diseases's interaction profile
for i in col0:
    for j in range(nm):
        Am_WNN[:,i] = Am_WNN[:,i] + (g**j)*A[:,argsort(-ms[i])[j]] #WNN infers miRNAs's interaction profile

A_score =zeros((ns,nm))
for i_M in range(M):
    SMS_subset=SMS[:,random.sample(range(ns),int(r*ns))]
    ms_subset=ms[:,random.sample(range(nm),int(r*nm))]
    SMS_subsetdimred = svdprocess(SMS_subset,int(r*ns))
    ms_subsetdimred = svdprocess(ms_subset,int(r*nm))
    
    SMS_subsetdimredkernel = zeros((SMS_subsetdimred.shape[0],SMS_subsetdimred.shape[0]))
    for i in range(SMS_subsetdimred.shape[0]):
        for j in range(SMS_subsetdimred.shape[0]):
            if j<=i:
                SMS_subsetdimredkernel[i,j]= exp (-(linalg.norm(SMS_subsetdimred[i]-SMS_subsetdimred[j])**2)/SMS_subsetdimred.shape[1])
    SMS_subsetdimredkernel = SMS_subsetdimredkernel + SMS_subsetdimredkernel.T-eye(ns)
    
    ms_subsetdimredkernel = zeros((ms_subsetdimred.shape[0],ms_subsetdimred.shape[0]))
    for i in range(ms_subsetdimred.shape[0]):
        for j in range(ms_subsetdimred.shape[0]):
            if j<=i:
               ms_subsetdimredkernel[i,j]= exp (-(linalg.norm(ms_subsetdimred[i]-ms_subsetdimred[j])**2)/ms_subsetdimred.shape[1])
    ms_subsetdimredkernel = ms_subsetdimredkernel + ms_subsetdimredkernel.T-eye(nm)        
    
    SM_regular=(mat(SMS_subsetdimredkernel)+regularpara*mat(eye(SMS_subsetdimredkernel.shape[0]))).I 
    m_regular=(mat(ms_subsetdimredkernel)+regularpara*mat(eye(ms_subsetdimredkernel.shape[0]))).I
    A_score= A_score + 0.5*mat(SMS_subsetdimredkernel)*SM_regular*mat(Am_WNN) + 0.5*(mat(ms_subsetdimredkernel)*m_regular*mat(As_WNN).T).T

A_scoremean = A_score/M
condition = A==0
A_0score = extract(condition,A_scoremean) #提取元素为0的得分
A_0scoreranknumber =  argsort(-A_0score)
dateset_n = argwhere(A == 0)
#SMrankname = []
#miRNArankname = []
for i in range(dateset_n.shape[0]):
    A_0scorerank = A_0score[A_0scoreranknumber[i]]
    SMrankname_pos = dateset_n[A_0scoreranknumber[i],0]
    SMrankname = SM_namenumber[SMrankname_pos,0]
    miRNArankname_pos = dateset_n[A_0scoreranknumber[i],1]
    miRNArankname = miRNA_namenumber[miRNArankname_pos,0]
    outscore = open(r'scorerank.csv', "a", newline = "")
    csv_writer = csv.writer(outscore, dialect = "excel")
    csv_writer.writerow([SMrankname,miRNArankname,A_0scorerank]) 
    outscore.close() 








       
      


