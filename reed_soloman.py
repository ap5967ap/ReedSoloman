import math
import random
import gmpy2
from helper import *

#*-------Global Variables-------------
M=cc('1'+'0'*1000)
mu=0
k=0
primes=[]

def GlobalSetup(mu,M):
    global k, primes
    k=cc(gmpy2.log10(M)+5)
    #print(k)
    primes=getK(k)
    primes.sort(reverse=1)

def ReedSolomonSend(a):
    global primes
    alist=[]
    for i in primes:
        alist.append(mod(a,i))
    return Transmit(alist)

def Transmit(a):
    global mu,k,primes
    l=rand(0,int(mu*k))
    I=set()
    while len(I)<l:
        temp=rand(0,k)
        if temp not in I:
            I.add(temp)
    # I=random.sample(range(k),random.randint(0,int(mu*k))) 
    b=[cc(0) for i in range(k)]
    for i in range(len(a)):
        if i in I:
            temp1=a[i]
            while temp1==a[i]:
                temp1=rand(cc(0),primes[i])
            b[i]=temp1
        else:
            b[i]=a[i]
    return b

def EEA(a,b):
    r=a
    rd=b
    s=cc(1)
    sd=cc(0)
    t=cc(0)
    td=cc(1)
    rlist=[r]
    slist=[s]
    tlist=[t]
    while(rd):
        q,rdd=r//rd, mod(r,rd)
        r,s,t,rd,sd,td=rd,sd,td,rdd,s-sd*q,t-td*q
        rlist.append(r)
        slist.append(s)
        tlist.append(t)
    return rlist,slist,tlist

def ReedSolomonReceive(b):
    global M,primes,mu,k
    P=cc(1)
    n=cc(1)
    l=int(mu*k)
    primes.sort(reverse=1) 
    for i in range(l):
        P=P*primes[i]
    for i in range(k):
        n=n*primes[i]
    # assert(n>2*M*P*P)
    bval=crt(primes,b)
    rlist,slist,tlist=EEA(n,bval)
    rst=M*P
    tst=P
    j=cc(0)
    for i in range(len(rlist)):
        if rlist[i]<=rst:
            j=i
            break
    rd=rlist[j]
    sd=slist[j]
    td=tlist[j]
    
    try:
        if (mod(rd,td))==0:
            return rd//td
        else:
            return -1
    except:
        return -1

def pre():
    global mu,M
    mu=float(input("Enter mu: \n=> "))
    M=cc(input("Enter the bound M: Enter 0 for default value of 10^1000\n=> "))
    if M==0:
        M=cc('1'+'0'*1000)
    GlobalSetup(mu,M)
    

   
def main():
    pre()
    t=int(input("Enter the number of test cases: \n=> "))
    for i in range(t):
        a=cc(input("Enter the message: \n=> "))
        b=ReedSolomonSend(a)
        print(f'Message {a} sent') #the input message
        print(f'Received message {ReedSolomonReceive(b)}') #the output message (-1 if unsuccessful reconstruction)

main()