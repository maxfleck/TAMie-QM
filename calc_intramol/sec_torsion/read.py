
import glob
import numpy as np
import matplotlib.pyplot as plt


hartree = 2625.5

files = glob.glob("tors_*/energy")

ens = []
deg = []
for f in files:
  f1 = open(f, "r")
  last_line = f1.readlines()[-2]
  num = float( [ x for x in last_line.split(" ") if x ][1] )
  ens.append(num)
  deg.append( float(f.split("/")[0].split("_")[1]) )
  f1.close()

ens = np.array(ens)*hartree
ens = ens - np.min(ens)

deg = np.array(deg)
p = np.argsort(deg)
deg = deg[p]
ens = ens[p]

plt.plot(deg,ens)
plt.savefig("torsion_energies.png")
plt.savefig("torsion_energies.pdf")
plt.show()
plt.close()


