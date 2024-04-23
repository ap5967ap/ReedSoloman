import random
import gmpy2
def power(a,b):
    c=gmpy2.mpz(1)
    n=a
    while(b>0):
        if b&1:
            c=c*a
        a=a*a
        b=b>>1
    return c

def powm(a,b,mod1):
    c=gmpy2.mpz(1)
    n=a
    while(b>0):
        if b&1:
            c=mod(c*a,mod1)
        a=mod(a*a,mod1)
        b=b>>1
    return c

def cc(a):
    return gmpy2.mpz(a)

def begcd(a,b):
    r=a
    rd=b
    e=gmpy2.mpz(0)
    while (r&1==0 and rd&1==0):
        r=r>>1
        rd=rd>>1
        e=e+1
    ad=r
    bd=rd
    s=gmpy2.mpz(1)
    t=gmpy2.mpz(0)
    sd=gmpy2.mpz(0)
    td=gmpy2.mpz(1)
    while rd:
        while r&1==0:
            r=r>>1
            if s&1==0 and t&1==0:
                s=s>>1
                t=t>>1
            else:
                s=(s+bd)>>1
                t=(t-ad)>>1
        while rd&1==0:
            rd=rd>>1
            if sd&1==0 and td&1==0:
                sd=sd>>1
                td=td>>1
            else:
                sd=(sd+bd)>>1
                td=(td-ad)>>1
        if r>rd:
            r,rd=rd,r
            s,sd=sd,s
            t,td=td,t
        rd=rd-r
        sd=sd-s
        td=td-t
    
    ans=cc(1)
    ans=ans<<e      
    ans=ans*r
    return ans,s,t

def modinv(a,mod):
    d,s,t=begcd(a,mod)
    assert(d==1)
    return s%mod

def mod(a,b):
    ans=a-b*(a//b)
    if ans<0:
        ans+=b
    return ans

def crt(primes_list,a_list):
    n=cc(1)
    for i in primes_list:
        n=n*i
    k=len(primes_list)
    nstar=[cc(1) for i in range(k)]
    e=[cc(1) for i in range(k)]
    for i in range(k):
        nstar[i]=n//primes_list[i]
        b=mod(nstar[i],primes_list[i])
        t=modinv(b,primes_list[i])
        e[i]=nstar[i]*t
    ans=cc(0)
    for i in range(k):
        ans=mod(ans+a_list[i]*e[i],n)
    return ans

def check(alpha,h,t,n):
    '''check if alpha is in the witness set for n'''
    beta=powm(alpha,t,n)
    if beta==cc(1):
        return True
    for i in range(h):
        if beta==n-1:
            return True
        if beta==cc(1):
            return False
        beta=mod(beta*beta,n)
    return False

def rand(a,b):
    x=min(a,b)
    y=max(a,b)
    a=x
    b=y
    return gmpy2.mpz(a+random.random()*(b-a))

def millar_rabin(n,threshold):
    if (n&1==cc(0) and n>cc(2)) or n==cc(1):
        return False
    if n==cc(2):
        return True
    nd=n-1
    ndp=n-1
    h=cc(0)
    while nd&1==0:
        nd=nd>>1
        h=h+1
    t=ndp//(cc(1)<<h)
    for i in range(threshold):
        a=rand(cc(2),n-1)
        if not check(a,h,t,n):
            return False
    return True

def getkprimes(k, lb, ub,threshold=10):
    primes=[]
    s=set()
    while len(primes)<k:
        a=rand(lb,ub)
        if millar_rabin(a,threshold) and a not in s:
            primes.append(a)
            s.add(a)
    return primes


def getK(k):
    primes=[]
    s=set()
    a=11
    while len(primes)<k:
        if millar_rabin(a,10) and a not in s:
            primes.append(a)
            s.add(a)
        a+=1
    return primes
    