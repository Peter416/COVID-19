import numpy as np
import matplotlib.pyplot as plt
N=14734014 # total population for normalizing the graphs
b=8.534e-9 # this beta seems to work reasonably well
a=0.0869 # same reason as b
n=150000 # vaccinated individuals per day
dt=0.001 # time scale between euler steps
t=np.arange(0,730,dt) # time vector for plotting
I=np.zeros((len(t))) # empty vector for storing values
S=np.zeros((len(t))) # same as above
R=np.zeros((len(t))) # same as above
cumul_I=np.zeros((len(t)))
cumul_R=np.zeros((len(t)))
S[0]=14353330 # initial conditions from April 1st 2020
I[0]=20875 #
R[0]=359809 #
print(np.size(I))
print(np.size(S))
print(np.size(R))
for x in range(1,len(t)): # for loop doing euler method of solving the System
  S[x]=dt*(-b*S[x-1]*I[x-1]-0.95*0.0007*S[x-1])+S[x-1] # euler step solutions
  I[x]=dt*(b*S[x-1]*I[x-1]-a*I[x-1])+I[x-1]
  R[x]=dt*(a*I[x-1]+0.95*0.0007*S[x-1])+R[x-1]
  if S[x]<0: # checks to see if susceptible goes negative
    S=[0]
    if I[x]<0.5: # checks to see if no infected are left
      I[x]=0 # If there are then the simulation ends
      S=S[:x]
      I=I[:x]
      R=R[:x]
      t=t[:x]
print('Pandemic has ended')
print(x)
break
plt.plot(t,S,'r',label='Susceptible') # all plots overlaid
plt.plot(t,I,'b',label='Infected')
plt.plot(t,R,'g',label='Removed')
plt.xlabel('time in days')
plt.ylabel('Population')
plt.title('SIR model for Onatrio starting from April 1st')
plt.legend()
plt.plot(S,I) #NOTE separate the blocks of code to produce graphs
plt.xlabel('S(t)')
plt.ylabel('I(t)')
plt.title('Phase plot for the SIR model')
