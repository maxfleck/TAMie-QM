from scipy.constants import epsilon_0,k,e
import numpy as np

def mie(r,p1,p2):
    # TO DO:
    # - add cut-off
    # - add shift

    # combining rules

    eps = np.sqrt(p1["epsilon"]*p2["epsilon"])
    
    if eps > 0:
        sig = (p1["sigma"]+p2["sigma"])/2
        m   = (p1["m"]+p2["m"])/2    
        # mie-force
        c0 = m / (m-6) * (m/6)**( 6/ (m-6) )
        mf = c0 * eps * ((sig / r) ** m - (sig / r) ** 6) 
        
    else:
        mf=0.0

    return mf

def mieq(r,p1,p2):
    # TO DO:
    # - add cut-off
    # - add shift
        
    # combining rules
    print("mie",r, p1["name"] , p2["name"] )
    eps = np.sqrt(p1["epsilon"]*p2["epsilon"])

    if eps > 0:
        sig = (p1["sigma"]+p2["sigma"])/2
        m   = (p1["m"]+p2["m"])/2    
        # mie-force
        c0 = m / (m-6) * (m/6)**( 6/ (m-6) )
        mf = c0 * eps * ((sig / r) ** m - (sig / r) ** 6) 
    else:
        mf=0.0
   
    # q-force
    #qf  = (p1["charge"]*p2["charge"] *1e10*e**2)/(4*np.pi*epsilon_0*r*k)
    if np.absolute(p1["charge"]) > 0.0 or np.absolute(p2["charge"] > 0):
        #print(e * e * 1.0e10 / ( 4.0 * np.pi * epsilon_0 * k ))
    
        qf = p1["charge"]*p2["charge"]*167100.002729011/r
    else:
        qf = 0.0
    return mf+qf


def qcharge(r,p1,p2):
    # TO DO:
    # - add cut-off
    # - add shift 
    # q-force
    print("ch",r, p1["name"] , p2["name"] )
    if np.absolute(p1["charge"]) > 0.0 or np.absolute(p2["charge"] > 0):
        qf = p1["charge"]*p2["charge"]*167100.002729011/r
    else:
        qf = 0.0    
    return qf


def bend(ang,p):
    # ang in radians!!
    if p["type"] == 2:
        force = 0.5*p["p"][1]*(ang-p["p"][0]*np.pi/180)**2  
    elif p["type"] == 3:
        force = 0.5*p["p"][1]*( np.cos(ang) - np.cos( p["p"][0]*np.pi/180 ) )**2
    else:
        force=0.0
    return force
    
def torsion(ang,pot):
    # ang in radians!!
    if isinstance(pot, dict):
        p = pot["p"]    
    elif isinstance(pot, list):
        p = pot
    elif isinstance(pot, np.ndarray):
        p = pot
    #print(ang)
    force = p[0] + p[1]* (1 + np.cos( ang )) + p[2]* (1 - np.cos(2* ang )) + p[3]* (1 + np.cos(3* ang ))
    
    return force
    
def torsion_types(ang,pot,key=1):
    # ang in radians!!
    if isinstance(pot, dict):
        p = pot["p"]    
    elif isinstance(pot, list):
        p = pot
    elif isinstance(pot, np.ndarray):
        p = pot
    #print(ang)
    if key == 1:
    	force = p[0] + p[1]* (1 + np.cos( ang )) + p[2]* (1 - np.cos(2* ang )) + p[3]* (1 + np.cos(3* ang ))
    elif key == 6:
    	force = p[0] + p[1]*np.cos( ang ) + p[2]*np.cos(ang)**2 + p[3]*np.cos(ang)**3 + p[4]*np.cos(ang)**4 + p[5]*np.cos(ang)**5 + p[6]*np.cos(ang)**6 + + p[7]*np.cos(ang)**7
    return force




