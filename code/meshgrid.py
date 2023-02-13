import numpy as np
from readData import Particle
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

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

    def __init__(self, x_min, x_max, y_min, y_max, square_count, Particle_list):
        self.x_min = x_min
        self.x_max = x_max
        self.y_min = y_min
        self.y_max = y_max
        self.length = x_max - x_min
        self.square_count = square_count
        self.h = self.length / self.square_count
        self.Particle_list = Particle_list


    def insert_particle(self):
        square_coor = {}
        part_list = []
        x_coor = np.arange(x_min,x_max+self.h,self.h)
        y_coor = np.arange(y_min,y_max+self.h,self.h)
        print(x_coor)
        print(y_coor)
        for i in range(len(x_coor)-1):

            for j in range(len(y_coor)-1):
                for particle in self.Particle_list:
                    if x_coor[i] <= particle.pos[0] <= x_coor[i+1]:
                        if y_coor[j] <= particle.pos[1] <= y_coor[j+1]:
                            part_list.append(particle.index)

                square_coor[round(x_coor[i],3),round(y_coor[j],3)] = part_list
                part_list = []

        return square_coor

    def plot(self):
        plt.scatter(particle_position_x, particle_position_y)
        # plt.scatter(particle_position_x, particle_position_y)
        plt.xlabel('x')
        plt.ylabel('y')
        plt.title('Particle Distribution')
        plt.show()



    def plot_histogram(self):
        particle_dict = self.insert_particle()

        particle_key_list = list(particle_dict.keys())
        print('YY')
        print(len(particle_key_list))
        data_matrix = np.array([len(particle_dict[i]) for i in particle_key_list])
        data_matrix = data_matrix.reshape(self.square_count,self.square_count)
        print(data_matrix)
        print(data_matrix.shape)


        plt.imshow(data_matrix, cmap='viridis', extent=[x_min,x_max,y_min,y_max])
        plt.title('Density Field')
        plt.colorbar()
        plt.show()









grid = Mesh(x_min,x_max,y_min,y_max,100,Particle_list)
grid.insert_particle()

grid.plot_histogram()





