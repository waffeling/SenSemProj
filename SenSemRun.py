#!/usr/bin/env python
# coding: utf-8

# In[14]:




##FOR E<VNAUGHT

import glob
import os
import datetime
import matplotlib.animation as animation
import matplotlib.pyplot as plt
from numpy import zeros, sqrt, sin, cos, pi, exp, cosh, sinh, format_float_scientific, linspace
from pylab import plot, legend
from matplotlib import rc
from IPython.display import HTML
##Trying again


def remove(String):
    String.replace(" ", "_")
    
#This is the total number of DISCRETE steps allowed in x 
Totalx = 1000


#These values should be adjusted should the chosen energies cause the program to blow up:
# Definition controls the size of an individual delx...
#   It's called definition as it should be adjusted to create appropriate definition in the initial wave
# Speed controls the size of an individual delt...
#   It's called speed as it should be adjusted to control the speed of the wave packet, and will also be the final
#value to adjust if the program continuously blows up (slower=better, sometimes)
definition = 0.04
speed = 0.0005






     

#Setting up our initial data and space
sigma = (Totalx*definition)/100 #Wavepacket size
start = (Totalx*definition)/4 #Wavepacket starting placement
eVmul = 0 #Ignore this
pause = False #Ignore this

#Some stupid global variables for the itegrator
Transm = 0
Refl = 0
integrate = 1
k = 0

#This is a cool function: it's what allows me to pause the animation
def onClick(event):
    global pause
    pause ^= True #makes pause
    global integrate
    integrate += 1 #this is how many times we want the integrator to work per pause (so it doesn't keep counting up)


    #THIS big block of code is meant to create all the units for the simulation
#It's based on prior math I did to force 1 kg' = 9.11e-31kg, hbar' = 1 J's' = 1.05e-34 Js, 
#and *hopefully* 1 J' = 1.6e-19 J
def unitcalculator(speed, definition, E, Vnaught):
    print("The following data is in reference to the current simulation: \n")

    # m' and s' are calculated first, as these are solely dependent on the prior unit considerations
    secondprime = (1.05*10**-34)/(1.6*10**(-19))
    specialnum = sqrt((1.6*10**(-19))/(9.11*10**-31))
    meterprime = secondprime * specialnum

    #k0 is then calculated under the previous unit considerations that E is in J' (eV) and k0 comes out in 1/m'
    #this should be true by hbar' = 1 J's' and kg' = 1me kg (I could use help checking this)
    k0 = sqrt(2*E)
    print("wave speed", k0, "mprime/sprime")

    #This data, as well as the above calculated data, is printed to ease the adjustment of definition and speed 
    #when a different magnitude of energies is to be tested
    wavelength = (2*pi)/k0
    print("wavelength", wavelength, "mprime")
    delx = definition
    print("delx", delx, "meterprime")
    delt = speed
    print("delt", delt, "secondprime")
    print("sim speed", delx/delt, "mprime/sprime")
    print("coeff", delt/(2*delx**2), "sprime/mprime^2")
    return meterprime, secondprime, k0, delx, delt, E, wavelength





#Setting up animation
fig = plt.figure()
ax = plt.axes(xlim=(0, Totalx-1), ylim=(0, 1))
ln1, = ax.plot([], [], lw=1)
ln2, = ax.plot([], [], lw=1, color='g')

#Creating labels to show calculated Transmitted and Reflected coefficients
time_template = 'T = %g'
time_template1 = 'R = %g'
time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)
time_text1 = ax.text(0.05, 0.85, '', transform=ax.transAxes)



#begin animation


def animate(j):
    global RealPsi
    global newRealPsi
    global V
    global newImPsi
    global ImPsi
    global FullPsi
    global Space
    global Time
    global Totalx
    global Totalt
    global blockbeg
    global blockend
    global k
    global Refl
    global fastfwd
    global pause
    global Transmlist
    global transmitted
    
    
     #allows the animation to be sped up for better viewing. Must be kept an EVEN INTEGER to allow leapfrogging
    if not pause: # allows pausing of the animation at any point for Transmition coefficicent calculation
        if fastfwd == 1:

            if (j%2)==0 and j != 0:

                newImPsi = zeros(Totalx)
                for i in range(Totalx):
                    if i == 0:
                        newImPsi[i] = 0
                    elif i == Totalx-1:
                        newImPsi[i] = 0
                    else:
                        newImPsi[i] = ImPsi[i] + (delt/(2*delx**2))*(RealPsi[i+1]-2*RealPsi[i]+RealPsi[i-1]) - delt*V[i]*RealPsi[i]
                ImPsi = newImPsi.copy()


            if (j%2) != 0 and j != 1:

                newRealPsi = zeros(Totalx)
                for i in range(Totalx):
                    if i == 0:
                        newRealPsi[i] = 0
                    elif i == Totalx-1:
                        newRealPsi[i] = 0
                    else:
                        newRealPsi[i] = RealPsi[i] - (delt/(2*delx**2))*(ImPsi[i+1]-2*ImPsi[i]+ImPsi[i-1]) + delt*V[i]*ImPsi[i]
                RealPsi = newRealPsi.copy()


            #putting together Psi^2
            for i in range(Totalx):
                Psi[i] = RealPsi[i]**2 + ImPsi[i]**2


        elif fastfwd > 1:
            for u in range(fastfwd):
                    if (u%2)==0:

                        newImPsi = zeros(Totalx)
                        for i in range(Totalx):
                            if i == 0:
                                newImPsi[i] = 0
                            elif i == Totalx-1:
                                newImPsi[i] = 0
                            else:
                                newImPsi[i] = ImPsi[i] + (delt/(2*delx**2))*(RealPsi[i+1]-2*RealPsi[i]+RealPsi[i-1]) - delt*V[i]*RealPsi[i]
                        ImPsi = newImPsi.copy()


                    if (u%2) != 0:

                        newRealPsi = zeros(Totalx)
                        for i in range(Totalx):
                            if i == 0:
                                newRealPsi[i] = 0
                            elif i == Totalx-1:
                                newRealPsi[i] = 0
                            else:
                                newRealPsi[i] = RealPsi[i] - (delt/(2*delx**2))*(ImPsi[i+1]-2*ImPsi[i]+ImPsi[i-1]) + delt*V[i]*ImPsi[i]
                        RealPsi = newRealPsi.copy()


                    #putting together Psi^2
                    for i in range(Totalx):
                        Psi[i] = RealPsi[i]**2 + ImPsi[i]**2

    else: 
        transmitted = 0
        reflected = 0
        total = 0
        if k < integrate:
            
            for i in range(Totalx):
                if i<blockbeg:
                    reflected += Psi[i]*delx
                if i>blockend:
                    transmitted += Psi[i]*delx
                if i<blockbeg or i>blockend:  
                    total += Psi[i]*delx
            k+=1
   

    
    if j >= (2*int(wavespinsim/delt)/fastfwd)-2 and not pause: 
        temptransmitted = 0
        reflected = 0
        total = 0
        for i in range(Totalx):
            if i<blockbeg:
                reflected += Psi[i]*delx
            if i>blockend:
                transmitted += Psi[i]*delx
            if i<blockbeg or i>blockend:  
                total += Psi[i]*delx
        
        Refl = reflected
        Transm = transmitted
        
        if transmitted < temptransmitted:
            transmitted = temptransmitted
            
            
        elif transmitted > temptransmitted:
            Transmlist.append(Transm)
            pause=True
            time_text.set_text(time_template % (Transm))
            time_text1.set_text(time_template1 % (Refl))
            
        

        
              

    
    
    ln2.set_data(Space, V)
    ln1.set_data(Space, Psi)
    
    

    return ln1, ln2, 


#this is what makes everything pause when you click on the animation
fig.canvas.mpl_connect('button_press_event', onClick)


#ignore these comments, they're for debugging purposes



#I know this bottom step makes the use of a function redundant, but it helps my brain organize the code
Totalx = 1000

E = 35
Eend = 65
Vnaught = 50
Transmlist = []
Elist = []
Estep = int(1)
Space = zeros(Totalx)
RealPsi = zeros(Totalx)
ImPsi = zeros(Totalx)
TotalPsi = 0
Psi = zeros(Totalx)
V = zeros(Totalx)
transmitted = 0

barrierlength = 2
a = barrierlength*(int(1/definition))

svdir = r"/home/pi/Desktop/" + str(datetime.date.today())
print(svdir)
remove(svdir)
print(svdir)

if os.path.exists(svdir) == False:
    os.makedirs(svdir)


while E < Eend:
    Primes = unitcalculator(speed, definition, E, Vnaught)
    meterprime = Primes[0]
    secondprime = Primes[1]
    k0 = Primes[2]
    delx = Primes[3]
    delt = Primes[4]
    E = Primes[5] 
    wavelength = Primes[6]

    extradis = 100
    blockbeg = 600
    blockend = blockbeg + a 
    wavespinsim = (blockend-start+extradis)*delx/k0
    #This is the total number of DISCRETE steps allowed in t 
    Totalt = 3*int(wavespinsim/delt)

    for i in range(Totalx):
        Space[i] = i
        RealPsi[i] = cos(k0*i*delx) * exp(-((i*delx)-start)**2/(sigma)**2) #turns out k0 must be multiplied by whatever delx you're using 
        ImPsi[i] = sin(k0*i*delx) * exp(-((i*delx)-start)**2/(sigma)**2)
        if i>= blockbeg and i<= blockend:
            V[i] = Vnaught
    TRealPsi = RealPsi.copy()
    TImPsi = ImPsi.copy()
    
    for i in range(Totalx):
        TotalPsi += (RealPsi[i]**2+ImPsi[i]**2)*delx
    
    for i in range(Totalx):
        RealPsi[i] = TRealPsi[i]/sqrt(TotalPsi)
        ImPsi[i] = TImPsi[i]/sqrt(TotalPsi)
    
    def init():
        ln1.set_data(Space, Psi)
        ln2.set_data(Space, V)
        return ln1, ln2, 

    pause = False
    fastfwd = 20

    anim = animation.FuncAnimation(
        fig, animate, init_func = init, interval=1, frames=int(Totalt/fastfwd), blit=True, repeat = True)
    f = svdir + r"/Eis" + str(E) + "Vis" + str(Vnaught) + ".mp4"
    writervideo = animation.FFMpegWriter(fps=60)
    anim.save(f, writer=writervideo)
    print(Transmlist)
    
    Elist.append(E)
    E += Estep

    
TrueData = zeros(100)
TestEs = linspace(0, Eend, 100)

u = 0

for i in TestEs:
    if i < Vnaught:
        TrueData[u] = (1/(1+((Vnaught**2)/(4*i*(Vnaught-i)))*sinh(sqrt(2*(Vnaught-i))*barrierlength)**2))
        u += 1
    elif i > Vnaught:
        TrueData[u] = (1/(1+((Vnaught**2)/(4*i*(i-Vnaught)))*sin(sqrt(2*(i-Vnaught))*barrierlength)**2))
        u += 1
    
print(TrueData)

TestEscp = TestEs.copy()
Elistcp = Elist.copy()
u=0

for i in TestEscp:
    TestEs[u]=TestEscp[u]/Vnaught
    u+=1
    
u=0
for i in Elist:
    Elist[u]=Elistcp[u]/Vnaught
    u+=1
    
anim.event_source.stop()

fig1 = plt.figure()
print(Elist)
plt.plot(Elist, Transmlist, "b-", label="Simulated Data")
plt.plot(TestEs, TrueData, "r-", label="Analytical Solution")
plt.xlim(E, 70)
plt.xlabel("E/V0")
plt.ylabel("T")
plt.legend(loc="upper left")
plt.title("Transmission Coefficient VS E/V0 | V0 = 50 and Barrier Width = 2")
plt.savefig(svdir + r'/Transmission_Coefficient_PlotEoverV.pdf')
