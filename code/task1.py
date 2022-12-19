import numpy as np
from readData import Particle
import math
import matplotlib.pyplot as plt

particles = np.genfromtxt("/Users/tunayildiz/Desktop/UZH/ComputationalAstrophysics/Data/data.txt", delimiter="\t", dtype=float) # ID, mass, x, y, z, vx, vy, vz, softening, potential
Particle_list = [Particle(i,m,x,y,z,vx,vy,vz,softening,potential) for i,m,x,y,z,vx,vy,vz,softening,potential in particles]
radius_max = int(max([i.radius for i in Particle_list])) + 1

def total_mass():
    mass = 0
    for i in Particle_list:
        mass += i.mass
    return mass


def Mass(r):  # Total mass in a given radius r.
    mass = 0
    mass_values = []
    for particle in Particle_list:
        if r >= particle.radius:
            mass += particle.mass

    return mass


def bin_mass(r, bin_size):
    mass = 0

    for particle in Particle_list:
        if particle.radius >= r and particle.radius <= r + bin_size:
            mass += particle.mass

    return mass

def volume(r,bin_size):

    return np.pi * pow(bin_size,2)


def half_mass(): # half mass radius = 0.188
    half = total_mass() / 2
    val = 0
    for i in np.arange(0, radius_max,0.001):
        if Mass(i) >= half:
            val = i
            break
    return val


scale_length = half_mass() / (1+math.sqrt(2))


def hernquist_density(r): # Hernquist Paper Equation 2.

    return (total_mass()/ (2*np.pi)) * (scale_length / r) * (1/pow((r+scale_length),3))


hernquist_vals = []
mass_vals = []
radius_list = []


for i in np.arange(0.01, radius_max,1):
    hernquist_vals.append(hernquist_density(i))
    mass_vals.append(bin_mass(i,5)/volume(i,5))


print(min(hernquist_vals))
print(max(hernquist_vals))


plt.plot(mass_vals,color = 'orange')
plt.show()




