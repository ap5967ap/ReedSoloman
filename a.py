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

M=cc('1'+'0'*17)
mu=random.random()
k=0
primes=[]

def GlobalSetup(mu,M):
    global k, primes
    k=16
    primes=getkprimes(k,cc(gmpy2.isqrt(M)),M-1)
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
        if rd%td==0:
            return rd//td
        else:
            return -1
    except:
        return -1
    
def main():
    GlobalSetup(mu,M)
    a=rand(0,M)
    b=ReedSolomonSend(a)
    print(a)
    # print(mu)
    print(ReedSolomonReceive(b))
    
if __name__ == '__main__':
    main()