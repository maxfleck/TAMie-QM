import numpy as np
import csv

class molecule():    
    def __init__(self,mol):
        self.molecule = np.array(mol)
        self.f = np.zeros(len(mol))
        self.n = -1*np.ones(len(mol))
        self.i = np.arange(len(mol))
        n=0
        for i,m in enumerate(mol):
            if m[0]=="b":
                self.f[i] = int(m[1:])
            elif m[0]=="r":
                self.f[i] = -int(m[1:])
            else:
                self.n[i] = n
                n+=1
        return
        
    def get_distance(self,nos,ignore_ring=False):
        n_min = np.min(nos)
        n_max = np.max(nos)
        i_min = np.squeeze(np.where(self.n==n_min))
        i_max = np.squeeze(np.where(self.n==n_max))
        #print(self.molecule[i_min],self.molecule[i_max])
        #print(self.molecule,i_min,i_max,n_min,n_max,nos)
        
        if i_min==0:
            seq=[]
            i_lastmol = i_min
            lastaction = 0
            flag=0
            #print(np.arange(i_min,i_max))
            for i in np.arange(i_min,i_max+1):
                if self.f[i]==0:
                    if flag==0:
                        seq.append(self.molecule[i])
                        i_lastmol = i
                        lastaction = 0
                    else:
                        flag-=1
                elif flag==0:
                    n_end = self.n[i_lastmol]+self.f[i]+lastaction
                    i_end = np.squeeze(np.where(self.n==n_end))

                    if i_end < i_max and self.f[i]>0:
                        flag = self.f[i] 
                        lastaction += self.f[i]
                    elif self.f[i]<0: #ring closing...here??
                        print("aaah")
                    #print(i,flag,i_end,i_max)
        else:
            #print(i_min,i_max,len(self.n))
            l1,s1=self.get_distance((0,int(n_min) ))   
            l2,s2=self.get_distance((0,int(n_max) ))   
            last=0
            for i,s in enumerate(s1):
                if s in s2:
                    last=i
                else:
                    break
            p1 = np.squeeze(np.where(s1==s1[last]))
            p2 = np.squeeze(np.where(s2==s2[last]))
            #print(last)
            #print("test",s1,s2)
            #print("test",p1,p2)
            #print(s1[p1:])
            #print(s2[p2:])
            ss1 = np.flip(s1[p1+1:])
            ss2 = s2[(p2):]
            seq= np.concatenate( (ss1,ss2) )
            
            
        return len(seq)-1,np.array(seq)
    
    def assign_xyz_coos(self,path):
        data = []
        with open(path) as csvfile:
            spamreader = csv.reader(csvfile,delimiter=" ")
            for row in spamreader:
                data.append([r for r in row if r])
        #print(data)
        n = int(data[0][0])
        xyz = {}
        for i,d in enumerate(data[2:n+2]):
            p = self.molecule[np.squeeze(np.where(self.n==i))]
            if d[0] in p.split("_")[0]:
                xyz[i] = {}
                xyz[i]["atom"] = p
                xyz[i]["original"] = p
                xyz[i]["xyz"] = np.array([float(x) for x in d[1:4]])    
            else:
                print("ERROR: xyz doesnt fit your structure")
                break
        return xyz


    def get_index(self,no):
        return np.squeeze(np.where(self.n==no))

