import numpy as np
import matplotlib.pyplot as plt 


p_infty = 101325.0
p_vap = 2300.0
p_ref = p_infty - p_vap
rho = 998.2067

r0 = 128e-6

t_ref = r0*np.sqrt(rho/p_ref)

Sigma = 0.07275/(r0*p_ref)
Pg0 = 1.0 + 2.0*Sigma
gamma = 7.0/5.0

print(t_ref)
print(Sigma)
print(Pg0)

def rp_force(x,v):
    return -1/x*(3/2* v**2 + 2*Sigma/x + 1 - Pg0*(1.0/x)**(3*gamma))

N = 100000

x = np.zeros(N)
v = np.zeros(N)

x[0] = 1.04
v[0] = 0.0

dt = 0.001

for i in range(1,N):
    v[i] = v[i-1] + rp_force(x[i-1],v[i-1])*dt
    x[i] = x[i-1] + v[i]*dt


t_rp = np.arange(N)*dt*t_ref*1000*1.0 # entspricht 15% decrease von Eigenfrequenz





data = np.fromfile('radius-time.txt',sep=';').reshape((-1,2))

fig = plt.figure()
ax = fig.add_subplot()

t = data[:,0]*1e3

x_f = np.sin(2*np.pi*30*t_rp)

x_i = x_f + (x-1)*25


ax.plot(t_rp,0.65 + 0.04*x_f)
ax.plot(t_rp,x-0.25)
ax.plot(t_rp,0.85 + 0.02*x_i)
ax.plot(t,data[:,1],'k-')
ax.set_xlabel('time [ms]')
ax.set_ylabel('radius []')

ax.text(-0.15,0.65,'p_acoustic')
ax.text(-0.15,0.75,'RP-sim')
ax.text(-0.15,0.85,'interference')
ax.text(-0.15,1.0,'3D-sim')

ax.set_xlim([-0.2,0.8])

plt.show()