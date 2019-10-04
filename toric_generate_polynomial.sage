 
# Takes a matrix g and a positive integer n and returns
# a list of powers of g, [g,g^2,...,g^n].
def gpows(g,n):
    out=[g]
    temp=g
    for i in range(n-1):
        temp=temp*g
        out.append(temp)
    return out
  
# Takes a sequence of k numbers (the traces of g acting on V)
# and returns the result of performing Newton-Girard on these nmbers
# (the traces of g acting on wedge 0,1,...,k of V).
def newtongirard(traces):
    trcs=[1]+traces
    k=len(trcs)
    trw=[0]*k
    trw[0]=trcs[0]
    for i in range(k-1):
        temp=0
        for j in range(i+1):
            temp=temp+(-1)**j*trcs[j+1]*trw[i-j]
        temp=temp/(i+1)
        trw[i+1]=temp
    return trw

#Computes the trace of the action of g on M.
def M_trace(g,M):
    out=0
    for i in range(M.matrix().dimensions()[0]):
        out=out+M.matrix().solve_left(g*M.matrix()[i])[i]
    return out

#Computes the number of torsion elements in Amb/M which are fixed by g.
def M_mult(g,M,Amb):
    Q=Amb/M
    invs=Q.invariants()
    i=0
    vecs=[vector(M.zero())]
    while i<len(invs) and invs[i]>1:
        v=Q.gen(i).lift()
        newvecs=[]
        for j in range(1,invs[i]):
            for k in range(len(vecs)):
                newvecs.append(vecs[k]+j*v)
        vecs=vecs+newvecs
        i=i+1
    out=0
    for v in vecs:
        if Q(v-g*v).is_zero():
            out=out+1
    return out

#Computes the value of the equivariant Poincaré polynomial of T=Hom(Amt/M,C^*) 
#at the group element g. Takes as inputs the relevant powers of g, their traces
#on the abient module and the module M.
def M_polynomial(gpws,Ambtraces,M):
    k=len(Ambtraces)-M.rank()
    Qtraces=[]
    for i in range(k):
       Qtraces.append(Ambtraces[i]-M_trace(gpws[i],M))
    wedges=newtongirard(Qtraces)
    t=var("t")
    pol=0
    for i in range(len(wedges)):
        pol=pol+wedges[i]*t^i
    return pol

#Computes the value of the equivariant Poincaré polynomial of the complement
#of the toric arrangement associated to the root system r at the group element g.
#Also takes (the underlying set of) the poset of intersections S and its Möbius
#function mu as input.
def torpol(g,r,S,mu):
    pol=0
    simp=r.ambient_space().simple_roots()
    Amb=[]
    for s in simp:
        Amb.append(vector(s))
    Amb=span(Amb,ZZ)
    d=Amb.rank()
    gpws=gpows(g,d)
    Ambtraces=[]
    t=var("t")
    for i in range(len(gpws)):
        Ambtraces.append(M_trace(gpws[i],Amb))
    for i in range(len(S)):
        if mod(i,25)==0:
            print("The polynomial is "+str((100*i/len(S)).n(digits=3))+"%"+" complete.")
        if mu[i]!=0:
            mult=M_mult(g,S[i],Amb)
            pol=pol+mult*mu[i]*M_polynomial(gpws,Ambtraces,S[i])*(-t)**(S[i].rank())
    return expand(pol)

#Computes values of the equivariant Poincaré polynomial at several elements.
def many_pols(els,r,start_ind,end_ind):
    out = []
    for i in range(start_ind,end_ind+1):
        temp = [i]
        s = stabsubmods(els[i],r)
        S = setofmodules(s)
        P = gen_rels(S)
        mu = combined_mob(P)
        pol = torpol(els[i],r,S,mu)
        temp.append(pol)
        out.append(temp)
    return out

# Computes values of the equivariant Poincaré polynomial for the action extended with {+/- 1}.
def pmpols(els,cls,mid,r,i):
    g = els[i]
    s = stabsubmods(g, r)
    S = setofmodules(s)
    P = gen_rels(S)
    mu = combined_mob(P)
    pol = torpol(g, r, S, mu)
    out = [[i,pol]]
    j = find_minus(g,mid,cls,0)
    if i !=j:
        mg = els[j]
        pol = torpol(mg, r, S, mu)
        out.append([j,pol])
    return out
