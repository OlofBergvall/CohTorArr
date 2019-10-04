# Converts the input of newmodules to something passable
# to parnewmodules.
def argtopar(S1,S2,S3):
    arg=[]
    for s in S1:
        arg.append((s,S2,S3))
    return arg
  
# Takes a module M (from the initial set of modules) and a set
# of previously generated modules (prev) and a total set of modules (old)
# and generates all sums of a module in prev and startM.
@parallel(sage.parallel.ncpus.ncpus())
def parnewmodules(M,prev,old):
    out=[]
    for P in prev:
        temp=M+P
        if temp not in old:
            temp=temp.echelonized_basis_matrix()
            temp.set_immutable()
            out.append(temp) #Apparantly we cannot parallelize using the type module as output, hence the weird conversion.
    return Set(out).list()

# Converts the output of parnewmodules to a set of modules.
def tomodlist(gen,pgend):
    out=[]
    for G in gen:
        for g in G[1]:
            #if isinstance(g,basestring):
            #    print(g)
            out.append(span(g,ZZ))
    return pgend.union(Set(out))

# Splits a list L into a bunch of lists with n elements.
# Since not all lists have a length divisible by L, the
# final list may be shorter than len(L)/n.
def split_list(L,n):
    l=(len(L)/n).floor()
    r=len(L)-n*l
    out=[]
    for i in range(0,l):
        out.append(Set([L[j] for j in range(i*n,(i+1)*n)]))
    if r>0:
        out.append(Set([L[j] for j in range(l*n,l*n+r)]))
    return out

# Takes an initial set of modules (start), a set
# of previously generated modules (prev) and a total set of modules (old)
# and generates all sums of a module in "prev" and a module in "start" which
# do not occur in "old".
def newmodules_par(start,prev,old):
    lists=split_list(start,Integer(sage.parallel.ncpus.ncpus()))
    nmds=Set([])
    for l in lists:
        arg=argtopar(l,prev,old)
        gen=parnewmodules(arg)
        nmds=tomodlist(gen,nmds)
    return nmds
  
# Takes a set of previously generated modules (prev), an initial 
# set of modules (start) and an total set of modules (old) and
# generates the next set of modules.
def newmodules(start,prev,old):
    out=Set([])
    for i in range(len(start)):
        for j in range(len(prev)):
            temp=start[i]+prev[j]
            if temp not in old:
                out=out.union(Set([temp]))
    return out
  
# Generates the set of g-stable submodules without poset
# structure.
def setofmodules(startlist):#, Integer n):
    t0=walltime()
    tcurr=walltime()
    start=startlist
    old=startlist
    prev=startlist
    ct=1
    while len(prev)>0:
        print(str(len(prev))+" sums of "+str(ct)+" modules were generated in "+str(walltime()-tcurr)+" seconds.")
        tcurr=walltime()
        prev=newmodules_par(start,prev,old)
        old=old+prev
        ct=ct+1
    n=old[0].degree()
    zv=vector([0]*n)
    ZM=old[0].span(zv)
    old=[ZM]+old.list()
    print("A set containing "+str(len(old))+" elements was generated in "+str(walltime()-t0)+" seconds.")
    return old
