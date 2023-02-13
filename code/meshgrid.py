import numpy as np
from readData import Particle
import matplotlib.pyplot as plt

particles = np.genfromtxt("/Users/tunayildiz/Desktop/UZH/ComputationalAstrophysics/Data/data1.txt", delimiter="\t", dtype=float)[1:] # ID, mass, x, y, z, vx, vy, vz, softening, potential

print(f'Number of particles: {len(particles)}')
Particle_list = [Particle(i,m,x,y,z,vx,vy,vz,softening,potential=0) for i,m,x,y,z,vx,vy,vz,softening in particles]


x_min = min([Particle.pos[0] for Particle in Particle_list])
x_max = max([Particle.pos[0] for Particle in Particle_list])
y_min = min([Particle.pos[1] for Particle in Particle_list])
y_max = max([Particle.pos[1] for Particle in Particle_list])
print(x_min,x_max,y_min,y_max)

particle_position_x = [i.pos[0] for i in Particle_list]
particle_position_y = [i.pos[1] for i in Particle_list]
particle_position_z = [i.pos[2] for i in Particle_list]

class Mesh:

    def __init__(self, x_min, x_max, y_min, y_max, square_count,Particle_list):
        self.x_min = x_min
        self.x_max = x_max
        self.y_min = y_min
        self.y_max = y_max
        self.length = x_max - x_min
        self.square_count = square_count
        self.h = self.length / self.square_count
        self.Particle_list = Particle_list

    def square_coordinates(self):
        square_coordinates = {}
        y_min = 0
        for idx_x, i in enumerate(np.arange(0,self.length,self.h)):
            x_min = 0

            for idx_y, j in enumerate(np.arange(0,self.length,self.h)):

                square_coordinates[idx_x,idx_y] = [x_min, x_min + self.h, y_min, y_min + self.h]
                x_min += self.h

            y_min += self.h

        return square_coordinates

    def insert_particle(self):
        particle_list = []
        particle_points = {}
        square_coordinates = self.square_coordinates()


        for key in square_coordinates:
            for particle in self.Particle_list:
                if square_coordinates[key][0] <= particle.pos[0] <= square_coordinates[key][1]:
                    if square_coordinates[key][2] <= particle.pos[1] <= square_coordinates[key][3]:
                        particle_list.append(particle.index)

            particle_points[key] = particle_list
            particle_list = []


        return particle_points

def plot():
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.scatter(particle_position_x,particle_position_y,particle_position_z)
    # plt.scatter(particle_position_x, particle_position_y)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Particle Distribution')
    plt.show()


grid = Mesh(x_min,x_max,y_min,y_max,50,Particle_list)
print(grid.insert_particle())







