import matplotlib.pyplot as plt
import numpy as np
'''
def plot_bl(data, filename):
  #x_range = (-0.002, 0.002)
  plt.hist(data[:,4], bins=50, color='skyblue', histtype='barstacked')#, range=x_range)
  #plt.xlim(x_range)
  plt.xlabel('z [m]', fontsize=12)
  plt.ylabel('particle numbers', fontsize=12)
  plt.title('number of electrons', fontsize=15)
  plt.savefig(filename, transparent='TURE')
  plt.close()

def plot_z_energy(data, filename):
    x_range = (-0.002, 0.002)
    plt.plot(data[:, 4], data[:, 5], linestyle='None', marker='o', markersize=6, alpha=0.2)
    plt.xlim(x_range)
    plt.xlabel('z [m]', fontsize=12)
    plt.ylabel('energy [$\gamma$]', fontsize=12)
    plt.title('Beam energy distribution', fontsize=15)
    plt.savefig(filename, transparent=True)
    plt.close()
'''
def plot_scan(data):
  plt.plot(data[:,0],data[:,1], '-bD', markersize=6, alpha=0.7)
  plt.xlabel('phase', fontsize=12)
  plt.ylabel('bunch length [m]', fontsize=12)
  #plt.ylabel('energy spread [$\gamma$]', fontsize=12)
  #plt.ylabel('energy [$\gamma$]', fontsize=12)
  plt.title('Scan Phase', fontsize=15)
  plt.savefig('scan phase', transparent='TURE')
  plt.close()
'''
def plot_energy(data,filename):
  plt.hist(data[:,5], bins=20, color='wheat')
  plt.xlabel('relative energy [MeV]', fontsize=12)
  plt.ylabel('particle number', fontsize=12)
  plt.title('Beam energy distribution', fontsize=15)
  plt.savefig(filename, transparent='TURE')
  plt.close()

def plot_xpx(data):
  plt.plot(data[:,0],data[:,1],'.', markersize=5 , label='x-px')
  plt.plot(data[:,2],data[:,3],'.', markersize=5, label='y-py')
  plt.legend(loc='best',fontsize=10,shadow=True,facecolor='#ccc',edgecolor='#000')
  plt.xlabel('x,y [m]', fontsize=12)
  plt.ylabel('px, py [mc]', fontsize=12)
  plt.title('phase space x, px, y, py', fontsize=15)
  plt.savefig('emittance', transparent='TURE')
  plt.close()

def plot_bl(a, b, c, d, e, f, g):
  
  plt.hist(a, bins=50, histtype='barstacked', alpha=0.7)
  plt.hist(b, bins=50, histtype='barstacked', alpha=0.7)
  plt.hist(c, bins=50, histtype='barstacked', alpha=0.7)
  plt.hist(d, bins=50, histtype='barstacked', alpha=0.7)
  plt.hist(e, bins=50, histtype='barstacked', alpha=0.7)
  plt.hist(f, bins=50, histtype='barstacked', alpha=0.7)
  plt.hist(g, bins=50, histtype='barstacked', alpha=0.7)
    
  plt.xlabel('z [m]', fontsize=12)
  plt.ylabel('particle numbers', fontsize=12)
  plt.title('bunch length', fontsize=15)
  plt.savefig('bunchlength', transparent='TURE')
  plt.close()
'''
scann = np.loadtxt('scanphase')

a = np.loadtxt('fort.20')
b = np.loadtxt('fort.21')
c = np.loadtxt('fort.22')
d = np.loadtxt('fort.23')
e = np.loadtxt('fort.24')
f = np.loadtxt('fort.25')
g = np.loadtxt('fort.26')

plot_scan(scann)

#plot_bl(a,'01')
#plot_bl(b,'02')
#plot_bl(c,'03')
#plot_bl(d,'04')
#plot_bl(e,'05')
#plot_bl(f,'06')
#plot_bl(g,'07')

#plot_energy(data2)
#plot_xpx(data3)