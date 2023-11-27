import numpy as np
import glob
import csv
from rapidfuzz.string_metric import normalized_levenshtein

# loose collection of functions doing somehing
# - read
# - write
# - angles, distances...
# - forces

def pair_from_row(row):
    pair={}
    pair["mass"] = float(row[1])
    pair["epsilon"] = float(row[2])
    pair["sigma"] = float(row[3])
    pair["m"] = float(row[4])
    pair["cut"] = float(row[5])
    pair["charge"] = 0.0
    return pair

def pair_of_h():
    pair={}
    pair["mass"] = 1.0
    pair["epsilon"] = 0.0
    pair["sigma"] = 1.0
    pair["m"] = 0
    pair["cut"] = 0.0
    pair["charge"] = 0.0
    return pair

def get_charge(row,pair_dict):
    keys = list(pair_dict.keys())
    dummy = np.zeros(len(keys))
    for i,key in enumerate(keys):
        dummy[i] = normalized_levenshtein(row[0][1:],key)
        #print(i, dummy[i])
    pointer = np.squeeze(np.where( dummy == np.amax(dummy)))
    charge = float(row[2])
    return keys[pointer],charge
    
def read_pair_potentials(path):
    keylist = ["!","#","models:","types","References"]
    pair_dict = {}
    flag=-1
    with open(path) as csvfile:
        spamreader = csv.reader(csvfile,delimiter=" ")
        for row in spamreader:
            row = [ x for x in row if x]
            if row:
                if "VdW-site" in row:
                    flag=1
                elif "coulomb-site" in row:
                    flag=2
                elif row[0] in keylist:
                    break            
                elif flag == 1 and row:
                    pair_dict[row[0]] = pair_from_row(row)
                elif flag == 2 and row:
                    #print(row)
                    if row[0].split("_")[0] == "cH":
                        pair_dict[row[0]] = pair_of_h()
                        pair_dict[row[0]]["charge"] = float(row[2])
                    else:
                        p,ch = get_charge(row,pair_dict) 
                        #print(p,ch)
                        pair_dict[p]["charge"] = ch
    return pair_dict
    
def read_bond_potentials(path):
    keylist = ["!","#","models:","types"]
    bond_dict = {}
    flag=0
    with open(path) as csvfile:
        spamreader = csv.reader(csvfile,delimiter=" ")
        for row in spamreader:
            row = [ x for x in row if x]
            if row:
                if row[0] == "model":
                    flag=1
                elif flag==1 and row[0] in keylist:
                    break
                elif flag==1:
                    name = "_".join(row[1:3])
                    bond_dict[name] = {}
                    bond_dict[name]["list"] = row[1:3]
                    bond_dict[name]["type"] = int(row[0]) 
                    #bond_dict[name]["len"] = float(row[3]) 
                    #bond_dict[name]["spring"] = float(row[4]) 
                    bond_dict[name]["p"] = [float(r) for r in row[3:]] 
                
    return bond_dict
    
def read_angle_potentials(path):
    keylist = ["!","#","models:","types"]
    angle_dict = {}
    flag=0
    with open(path) as csvfile:
        spamreader = csv.reader(csvfile,delimiter=" ")
        for row in spamreader:
            row = [ x for x in row if x]
            if row:
                if row[0] == "model":
                    flag=1
                elif flag==1 and row[0] in keylist:
                    break
                elif flag==1:
                    name = "_".join(row[1:4])
                    angle_dict[name] = {}
                    angle_dict[name]["list"] = row[1:4]
                    angle_dict[name]["type"] = int(row[0]) 
                    #angle_dict[name]["angle"] = float(row[4]) 
                    #angle_dict[name]["p"] = float(row[5]) 
                    angle_dict[name]["p"] = [float(r) for r in row[4:]] 
                
    return angle_dict

def read_torsion_potentials(path):
    keylist = ["!","#","models:","types","Note:"]
    torsion_dict = {}
    flag=0
    with open(path) as csvfile:
        spamreader = csv.reader(csvfile,delimiter=" ")
        for row in spamreader:
            row = [ x for x in row if x]
            if row:
                if row[0] == "model":
                    flag=1
                elif flag==1 and row[0] in keylist:
                    break
                elif flag==1:
                    name = "_".join(row[1:5])
                    torsion_dict[name] = {}
                    torsion_dict[name]["list"] = row[1:5]
                    torsion_dict[name]["type"] = int(row[0]) 
                    torsion_dict[name]["p"] = [float(r) for r in row[5:]] 
                
    return torsion_dict

    
def read_xyz(path,energy=False):
    data = []
    with open(path) as csvfile:
        spamreader = csv.reader(csvfile,delimiter=" ")
        for row in spamreader:
            data.append([r for r in row if r])
    n = int(data[0][0])
    xyz = {}
    for i,d in enumerate(data[2:n+2]):
        xyz[i] = {}
        xyz[i]["atom"] = d[0]
        xyz[i]["xyz"] = np.array([float(x) for x in d[1:4]])          

    if energy:
        energy = float(data[1][2])
        return xyz,energy
        
    return xyz

def assign_CHx(xyz):
    xyz_CHx = {}
    xkeys = sorted(xyz.keys())
    ii=0
    x=np.min(xkeys)
    flag=0
    while x <= np.max(xkeys):
        if flag==0 and xyz[x]["atom"] == "C":
            coos = np.array(xyz[x]["xyz"])
            flag=1
            x+=1
        elif flag>0:
            if xyz[x]["atom"] == "H":
                coos = np.column_stack((coos, xyz[x]["xyz"]))
                x+=1
                flag+=1
            else:
                no=flag-1
                xyz_CHx[ii] = {}
                xyz_CHx[ii]["atom"] = "CH"+str(no)
                xyz_CHx[ii]["xyz"] = np.sum(coos*np.array([15]+[1]*no),axis=1)/(15+1*no)
                #print(xyz_CHx[ii]["xyz"])
                coos = []
                flag=0
                ii+=1
        else:
            xyz_CHx[ii] = xyz[x]
            x+=1
            ii+=1
    return xyz_CHx

def kill_CHx(xyz):
    xyz_CHx = {}
    xkeys = sorted(xyz.keys())
    ii=0
    x=np.min(xkeys)
    flag=0
    while x <= np.max(xkeys):
        if flag==0 and xyz[x]["atom"] == "C":
            coos = np.array(xyz[x]["xyz"])
            flag=1
            x+=1
        elif flag>0:
            if xyz[x]["atom"] == "H":
                #coos = np.column_stack((coos, xyz[x]["xyz"]))
                x+=1
                flag+=1
            else:
                no=flag-1
                xyz_CHx[ii] = {}
                xyz_CHx[ii]["atom"] = "C"
                xyz_CHx[ii]["xyz"] = coos
                #print(xyz_CHx[ii]["xyz"])
                coos = []
                flag=0
                ii+=1
        else:
            xyz_CHx[ii] = xyz[x]
            x+=1
            ii+=1
    return xyz_CHx


        
def to_xyz(xyz,path):
    f = open(path, "w")
    f.write(str(len(xyz))+"\n")
    f.write("\n")
    for x in xyz:
        line = "    ".join([xyz[x]["atom"][0]]+[str(y) for y in xyz[x]["xyz"]])
        f.write(line+"\n")
    f.close()
    
    
    
def distance(x1,x2):
    # distances
    r = np.linalg.norm(x1-x2)
    #print("r:",r)
    return r

def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)

def angle_between(x1, x2, x3):
    """ Returns the angle in radians between vectors 'v1' and 'v2'::

            >>> angle_between((1, 0, 0), (0, 1, 0))
            1.5707963267948966
            >>> angle_between((1, 0, 0), (1, 0, 0))
            0.0
            >>> angle_between((1, 0, 0), (-1, 0, 0))
            3.141592653589793
    """
    v1 = unit_vector(x1-x2)
    v2 = unit_vector(x3-x2)
    
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

def dihedral(x0,x1,x2,x3):
    """Praxeolitic formula
    1 sqrt, 1 cross product"""

    b0 = -1.0*(x1 - x0)
    b1 = x2 - x1
    b2 = x3 - x2

    # normalize b1 so that it does not influence magnitude of vector
    # rejections that come next
    #print(b1)
    #print(np.linalg.norm(b1))
    if np.linalg.norm(b1)>0:
        b1 /= np.linalg.norm(b1)

    # vector rejections
    # v = projection of b0 onto plane perpendicular to b1
    #   = b0 minus component that aligns with b1
    # w = projection of b2 onto plane perpendicular to b1
    #   = b2 minus component that aligns with b1
    v = b0 - np.dot(b0, b1)*b1
    w = b2 - np.dot(b2, b1)*b1

    # angle between v and w in a plane is the torsion angle
    # v and w may not be normalized but that's fine since tan is y/x
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return np.degrees(np.arctan2(y, x))
    
def read_json(path):
    with open(path) as json_file:
        return json.load(json_file)
