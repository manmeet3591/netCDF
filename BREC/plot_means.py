mport pylab as pl
from pylab import genfromtxt

cya = genfromtxt("cya.txt")
lpp = genfromtxt("lpp.txt")
spp = genfromtxt("spp.txt")
zoo = genfromtxt("zoo.txt")
nh4 = genfromtxt("nh4.txt")
salt = genfromtxt("salt.txt")
temp = genfromtxt("temp.txt")
det = genfromtxt("det.txt")

pl.figure(1)
pl.plot(cya[:,0], cya[:,1])
pl.title('Cya over years')

pl.figure(2)
pl.plot(temp[:,0], temp[:,1])
pl.title('Temp over years')

pl.figure(3)
pl.plot(lpp[:,0], lpp[:,1])
pl.title('lpp over years')

pl.figure(4)
pl.plot(spp[:,0], spp[:,1])
pl.title('spp over years')

pl.show()
