from a import *
M=cc('1'+'0'*200)
l=[]
i=gmpy2.isqrt(M)
while len(l)<15:
    if millar_rabin(i,10) and i not in l:
        l.append(i)
    i+=1
print(l)