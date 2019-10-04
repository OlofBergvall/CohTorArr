# Converts the input of gen_rels to something passable
# to par_gen_rels.
def argtopar_gen_rels(L,S):
    arg=[]
    for l in L:
        arg.append((l,S))
    return arg

# Computes inclusion relations of the poset of modules.
@parallel(sage.parallel.ncpus.ncpus())
def par_gen_rels(M,S):
    m=matrix(1,len(S),sparse=True)
    #m=[0 for s in S]
    for i in range(len(S)):
        if M.is_submodule(S[i]):
            #m[i]=1
            m[0,i]=1
    return m

# Converts the output of parnewmodules to a set of modules.
def tolists(gen,pgendrel,pgendel):
    outrel=[]
    outel=[]
    for G in gen:
        outrel.append(G[1])
        outel.append(G[0][0][0])
    outrel=block_matrix(len(outrel),1,outrel)
    if pgendrel!=[]:
        outrel=block_matrix(2,1,[pgendrel,outrel])
    outel=pgendel+outel
    return [outrel,outel]

# Splits a list L into a bunch of lists with n elements.
# Since not all lists have a length divisible by L, the
# final list may be shorter than len(L)/n.
def split_list(L,n):
    l=(len(L)/n).floor()
    r=len(L)-n*l
    out=[]
    for i in range(0,l):
        out.append([L[j] for j in range(i*n,(i+1)*n)])
    if r>0:
        out.append([L[j] for j in range(l*n,l*n+r)])
    return out

# Takes two lists L and S with the same entries but in different order
# and returns a list of integers such that the i'th entry is the index
# of L[i] in S.
def give_inds(L,S):
    out=[]
    for l in L:
        b=l==S[0]
        i=0
        while b==False:
           i=i+1
           b=l==S[i]
        out.append(i)
    return out

# Takes a list L and a list of integers inds and returns
# a list with the same entries as L but with the order
# specified by inds.
def perm_list(L,inds):
    out=[0 for i in inds]
    for i in range(len(inds)):
        out[inds[i]]=L[i]
    return matrix(out,sparse=True)
            
# Takes a set S of submodules of an ambient module M 
# and returns the inclusion relations among the elements of S.
def gen_rels(S):
    t0=walltime()
    rels=[]
    els=[]
    L=split_list(S,Integer(10*sage.parallel.ncpus.ncpus()))
    ct=0
    for l in L:
        print("The relations of the poset have been "+str((100*ct/len(L)).n(digits=3))+"%"+" completed in "+str(walltime()-t0)+" seconds.")
        arg=argtopar_gen_rels(l,S)
        gen=par_gen_rels(arg)
        temp=tolists(gen,rels,els)
        rels=temp[0]
        els=temp[1]
        ct=ct+1
    inds=give_inds(els,S)
    rels=perm_list(rels,inds)
    print("The relations of the poset were computed in "+str(walltime()-t0)+" seconds.")
    return rels

# Help function for printing date.
def add_zero(n):
    if n<10:
        return str(0)+str(n)
    else:
        return str(n)

#Returns string with date information.
def give_date():
    from datetime import datetime
    now=datetime.now()
    return str(now.year)+"-"+add_zero(now.month)+"-"+add_zero(now.day)+", at "+add_zero(now.hour)+":"+add_zero(now.minute)+":"+add_zero(now.second)
    
#Generates a subset of the inclusion relations (in case one wants to
#parallelize over several computers/nodes).
def part_gen_rels(S,stind,endind,list_name,file_name):
    Spart=[S[i] for i in range(stind,endind+1)]
    t0=walltime()
    rels=[]
    els=[]
    L=split_list(Spart,Integer(10*sage.parallel.ncpus.ncpus()))
    ct=0
    for l in L:
        print(give_date())
        print("The relations of the poset has been "+str((100*ct/len(L)).n(digits=3))+"%"+" completed in "+str(walltime()-t0)+" seconds.")
        arg=argtopar_gen_rels(l,S)
        gen=par_gen_rels(arg)
        temp=tolists(gen,rels,els)
        rels=temp[0]
        els=temp[1]
        ct=ct+1
    inds=give_inds(els,Spart)
    rels=perm_list(rels,inds)
    rels=sparse_zero_one_out(rels)
    print("A list of "+str(endind-stind+1)+" relations was computed in "+str(walltime()-t0)+" seconds.")
    save_rels(rels,list_name,file_name)
    return rels

#Repeats part_gen_rels N times (to avoid having to restart computations
# but keeping sizes of files reasonably big).
def many_gen_rels(S,stind,endind,N,Rstr):
    n=(endind-stind+1)/N
    for i in range(1,N+1):
        si=stind+(i-1)*n
        ei=stind+i*n-1
        endstr=str(si)+"_"+str(ei)
        ln="R_"+endstr
        fn=Rstr+"_rels_"+endstr
        rels=part_gen_rels(S,si,ei,ln,fn)
    print("Done!")

#Three methods for computing Möbius functions of a poset.
#The first two takes a list of relations as input, the third
#expects the poset to be given in the form of a matrix.
def mob(rel_list):
    M=matrix(rel_list)
    v=matrix(1,len(rel_list))
    v[0,0]=1
    return M.solve_left(v).list()
  
def mob2(rel_list):
    M=block_matrix(len(rel_list),rel_list)
    v=matrix(1,len(rel_list))
    v[0,0]=1
    return M.solve_left(v).list()
  
def combined_mob(P):
    v=matrix(1,P.dimensions()[0])
    v[0,0]=1
    return P.solve_left(v).list()

# Method for saving relations in a longer computation.
def save_rels(rels,list_name,file_name):
    print("Saving relations to \""+list_name+"\" in the file \""+file_name+".sage\".")
    f=open(file_name+".sage","w")
    f.write(list_name+"="+str(rels))
    f.close()
    print("Saving complete.")
    
# Converts a matrix of zeros and ones to the form
# [nrows,ncols,[list of column indicies where row 1 ha a 1],...,[list of column indicies where row nrows ha a 1]]
def sparse_zero_one_out(M):
    m=M.dimensions()[0]
    n=M.dimensions()[1]
    out=[m,n]
    for i in range(m):
        temp=[]
        for j in range(n):
            if M[i,j]==1:
                temp.append(j)
        out.append(temp)
    return out

# Converts the output of part_gen_rels into a zero-one matrix
# (i.e. converts inclusion information given in the form "element i is contained in elements [j1,j2,...,jm]"
# to inclusion information provided by the corresponding zero-one matrix).
def sparse_zero_one_in(L):
    m=L[0]
    n=L[1]
    out=matrix(m,n)
    for i in range(n):
        r=L[i+2]
        for j in range(len(r)):
            out[i,r[j]]=1
    return out
