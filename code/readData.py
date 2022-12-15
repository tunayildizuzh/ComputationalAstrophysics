import numpy as np


class Particle():

    def __init__(self,index,mass,x,y,z,vx,vy,vz,softening,potential):
        self.index = index
        self.mass = mass
        self.x = x
        self.y = y
        self.z = z
        self.vx = vx
        self.vy = vy
        self.vz = vz
        self.softening = softening
        self.potential = potential
        self.pos = np.array(self.x, self.y, self.z)
        self.radius = self.pos / np.linalg.norm(self.pos)


particles = np.genfromtxt("/Users/tunayildiz/Desktop/UZH/ComputationalAstrophysics/Data/data.txt", delimiter="\t", dtype=float) # ID, mass, x, y, z, vx, vy, vz, softening, potential

Particle_list = [Particle(i,m,x,y,z,vx,vy,vz,softening,potential) for i,m,x,y,z,vx,vy,vz,softening,potential in particles]

def total_mass():
    mass = 0
    for i in Particle_list:
        mass += i.mass
    return mass


