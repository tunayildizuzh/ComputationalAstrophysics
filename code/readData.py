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
        self.pos = np.array((self.x, self.y, self.z))
        self.radius = np.linalg.norm(self.pos)
        self.radius2 = np.linalg.norm(self.pos, ord = 2) ** 2
        self.force = 0




