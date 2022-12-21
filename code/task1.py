import numpy as np
from readData import Particle
import math
import matplotlib.pyplot as plt

particles = np.genfromtxt("/Users/tunayildiz/Desktop/UZH/ComputationalAstrophysics/Data/data.txt", delimiter="\t", dtype=float) # ID, mass, x, y, z, vx, vy, vz, softening, potential
Particle_list = [Particle(i,m,x,y,z,vx,vy,vz,softening,potential) for i,m,x,y,z,vx,vy,vz,softening,potential in particles]
radius_max = int(max([i.radius for i in Particle_list])) + 1

force_100 = np.genfromtxt('/Users/tunayildiz/Desktop/UZH/ComputationalAstrophysics/Data/nbody_forces100.txt', delimiter ='\t', dtype = float)
force_half = np.genfromtxt('/Users/tunayildiz/Desktop/UZH/ComputationalAstrophysics/Data/nbody_forces0.5.txt',delimiter = '\t', dtype = float)
force_rhm = np.hsplit(np.genfromtxt('/Users/tunayildiz/Desktop/UZH/ComputationalAstrophysics/Data/nbody_forces_rhm.txt',delimiter = '\t', dtype = float),2)[1]
particle_radius = [i.radius for i in Particle_list]

def total_mass():
    mass = 0
    for i in Particle_list:
        mass += i.mass
    return mass


def Mass(r):  # Total mass in a given radius r.
    mass = 0
    count = 0
    for particle in Particle_list:
        if r >= particle.radius:
            mass += particle.mass
            count += 1

    return mass, count


def bin_mass(r, bin_size):
    mass = 0

    for particle in Particle_list:
        if particle.radius >= r and particle.radius <= r + bin_size:
            mass += particle.mass

    count = mass / Particle_list[0].mass # Every particle has the same mass.
    return mass, count


def volume(r,bin_size):

    return np.pi * pow(bin_size,2)


def half_mass(): # half mass radius = 0.188
    half = total_mass() / 2
    val = 0
    for i in np.arange(0, radius_max,0.001):
        if Mass(i)[0] >= half:
            val = i
            break
    return val


softening = Mass(half_mass()/2)[1]
scale_length = half_mass() / (1+math.sqrt(2))


def hernquist_density(r): # Hernquist Paper Equation 2.

    return (total_mass()/ (2*np.pi)) * (scale_length / r) * (1/pow((r+scale_length),3))


hernquist_vals = []
mass_vals = []
radius_list = []
particle_count = []

for i in np.arange(0.01, radius_max,1):
    hernquist_vals.append(hernquist_density(i))
    mass_vals.append(bin_mass(i,5)[0]/volume(i,5))
    particle_count.append(bin_mass(i, 5)[1])

def compute_nbody_forces(G, epsilon):
    print('Computing N-body Forces')
    forces = []
    i = 0
    with open('nbody_forces0.5.txt', 'w') as output:  # Save it in text file.
        for particle in Particle_list:
            i += 1
            if i % 1000 == 0:
                print(i)
            r = 0
            for particle2 in Particle_list:
                r += particle2.mass / pow((pow(particle.radius - particle2.radius, 2) + pow(epsilon, 3)), (3 / 2)) * (
                            particle.radius - particle2.radius)

            output.write(str(-G*r) + '\n')
        output.close()

    return np.array(forces)

def analytical_force(G):
    force_analytic = []
    for i in np.arange(0.01, radius_max/350,1):
        force_analytic.append(((-G) * Mass(i)[0]) / pow(i, 3))

    return force_analytic


def draw_figs(number):  # Plots
    if number == 1:
        fig, (ax1,ax2,ax3) = plt.subplots(1,3, figsize =(15,6))
        ax1.plot(mass_vals)
        ax1.set_title('Mass Density Function')
        ax1.set_xlabel('Radius')
        ax1.set_ylabel('Mass Density')

        ax2.plot(hernquist_vals)
        ax2.set_title('Analytical Hernquist Density Function')
        ax2.set_xlabel('Radius')
        ax2.set_ylabel('Hernquist Density')

        ax3.plot(particle_count)
        ax3.set_title('Particle Count')
        ax3.set_xlabel('Radius')
        ax3.set_ylabel('Particle Count')
        plt.savefig('HernquistvsMassDensityFunction.png')
        plt.show()

    if number == 2:
        fig, (ax1,ax2,ax3) = plt.subplots(1,3, figsize = (12,8))
        ax1.plot(particle_radius,force_100, label = 'softening = 100', color = 'blue')
        ax2.plot(force_rhm, label='softening = rhm', color='orange')
        ax3.plot(particle_radius,force_half, label = 'softening = rhm/2', color = 'green')

        for i in ax1,ax3:
            i.set_xlabel('Particle Index')
            i.set_ylabel('Force')
        ax1.legend()
        ax2.legend()
        ax3.legend()
        plt.savefig('ForcevsSoftening.png')
        plt.show()

    if number == 3:
        plt.plot(analytical_force(1))
        plt.show()






# compute_nbody_forces(1, softening)
draw_figs(3)





