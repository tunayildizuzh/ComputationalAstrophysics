import numpy as np
from readData import Particle
import math
import matplotlib.pyplot as plt

particles = np.genfromtxt("/Users/tunayildiz/Desktop/UZH/ComputationalAstrophysics/Data/data_neww.txt", delimiter="\t", dtype=float) # ID, mass, x, y, z, vx, vy, vz, softening, potential
print(f'Number of particles: {len(particles)}')
Particle_list = [Particle(i,m,x,y,z,vx,vy,vz,softening,potential) for i,m,x,y,z,vx,vy,vz,softening,potential in particles]
radius_max = int(max([i.radius for i in Particle_list])) + 1


particle_radius = [i.radius for i in Particle_list]
particle_position_x = [i.pos[0] for i in Particle_list]
particle_position_y = [i.pos[1] for i in Particle_list]
particle_position_z = [i.pos[2] for i in Particle_list]

def total_mass():
    mass = 0
    for i in Particle_list:
        mass += i.mass
    return mass


def Mass(r):  # Total mass in a given radius r and number of particles in r.
    mass = 0
    count = 0
    for particle in Particle_list:
        if r >= particle.radius:
            mass += particle.mass
            count += 1

    return mass, count


def bin_mass(r, bin_size): # Mass and number of particals in a given shell.
    mass = 0

    for particle in Particle_list:
        if particle.radius >= r and particle.radius <= r + bin_size:
            mass += particle.mass

    count = mass / Particle_list[0].mass # Every particle has the same mass.
    return mass, count


def volume(r,bin_size):

    return np.pi * pow(bin_size,2)


def half_mass_radius(): # half mass radius = 0.188
    half = total_mass() / 2
    val = 0
    for i in np.arange(0, radius_max,0.001):
        if Mass(i)[0] >= half:
            val = i
            break
    return val


softening = Mass(half_mass_radius()/2)[1]
softening_mean_inter_sep = (pow(half_mass_radius(), 3) / pow(len(particles), (1/3)))
print(f'Mean Interparticle Seperation: {softening_mean_inter_sep}')
scale_length = half_mass_radius() / (1+math.sqrt(2))


def hernquist_density(r): # Hernquist Paper Equation 2.

    return (total_mass()/ (2*np.pi)) * (scale_length / r) * (1/pow((r+scale_length),3))


def compute_density_comparison():
    hernquist_vals = []
    mass_density_vals = []
    particle_count = []

    print('----- TASK1 -----')
    print('Hernquist Density and Mass Density calculation starts.')
    for i in np.arange(0.01, radius_max,1):
        hernquist_vals.append(hernquist_density(i))
        mass_density_vals.append(bin_mass(i,5)[0]/volume(i,5))
        particle_count.append(bin_mass(i, 5)[1])
    print('Hernquist Density and Mass Density calculation ended.')

    return hernquist_vals,mass_density_vals,particle_count


# densities = compute_density_comparison() # Hernquist, Mass, particle count



def nbody_force(particlei, particlej, epsilon):
    sum_r = 0
    if particlei.index != particlej.index:
        diff = particlei.pos - particlej.pos
        d2 = np.sum(diff**2)
        sum_r += particlej.mass / pow((d2 + epsilon**2), (3/2)) * diff
    force = (-1) * sum_r * particlei.mass
    force_magnitude = np.linalg.norm(force)
    return force_magnitude


def compute_nbody_forces2(G, epsilon,fname):
    print('----- TASK2 -----')
    print('Computing N-body Forces')
    forces = []
    i = 0
    with open(fname, 'w') as output:  # Save it in text file.
        for particlei in Particle_list:
            i += 1
            if i % 1000 == 0:
                print(f'{i} out of {len(particles)} calculated.')
            for particlej in Particle_list:
                force_magnitude = nbody_force(particlei, particlej, epsilon)


            particlei.force = force_magnitude
            output.write(str(-G*force_magnitude) + '\n')
        output.close()

# compute_nbody_forces2(1,softening_mean_inter_sep, 'direct_nbody_forces_mis_new.txt')
# compute_nbody_forces2(1,softening_mean_inter_sep*200,'direct_nbody_forces_mis2x_new.txt')
# compute_nbody_forces2(1,softening_mean_inter_sep*400,'direct_nbody_forces_mis0.5x_new.txt')
#
# force_mis_new = np.genfromtxt('/Users/tunayildiz/Desktop/UZH/ComputationalAstrophysics/code/direct_nbody_forces_mis_new.txt',delimiter ='\t', dtype=float)
# force_mis_new2 = np.genfromtxt('/Users/tunayildiz/Desktop/UZH/ComputationalAstrophysics/code/direct_nbody_forces_mis2x_new.txt',delimiter ='\t', dtype=float)
# force_mis_new3 = np.genfromtxt('/Users/tunayildiz/Desktop/UZH/ComputationalAstrophysics/code/direct_nbody_forces_mis0.5x_new.txt',delimiter ='\t', dtype=float)



# def compute_nbody_forces(G, epsilon):
#     print('----- TASK2 -----')
#     print('Computing N-body Forces')
#     forces = []
#     i = 0
#     with open('direct_nbody_forces_mis.txt', 'w') as output:  # Save it in text file.
#         for particle in Particle_list:
#             i += 1
#             if i % 1000 == 0:
#                 print(i)
#             r = 0
#             for particle2 in Particle_list:
#                 r += particle2.mass / pow((pow(particle.radius - particle2.radius, 2) + pow(epsilon, 2)), (3 / 2)) * (
#                             particle.radius - particle2.radius)
#
#             output.write(str(-G*r) + '\n')
#         output.close()
#
#         for idx, val in enumerate(Particle_list):
#             val.force = forces[idx]
#
#     return np.array(forces)


def analytic_force(G,r):
    return -((G*Mass(r)[0])/pow(r,3))


def compute_analytical_force():
    analytical_force_list = []
    print('Analytical Force Calculation starts.')

    for i in np.arange(0.01, radius_max,1):
        analytical_force_list.append(analytic_force(-1,i))

    print('Analytical Force Calculation ended.')
    return analytical_force_list





def compute_relaxation():
    G = 1
    N = len(Particle_list)
    vc = math.sqrt((G * total_mass()/2) / half_mass_radius()) #  M(Rhm) = total_mass() / 2
    t_cross = half_mass_radius() / vc
    t_relax = (N / (8 * np.log(N))) * t_cross
    # print(G,N,vc,t_cross,t_relax) # 1 50010 3496.8690678643025 5.404834906084457e-05 0.031226471422577867
    return G,N,vc,t_cross,t_relax

def leap_frog(run_time, dt):
    t = 0
    while t < run_time:
        for particle in Particle_list:
            particle.velocity += particle.force * dt/2.0 # 1/2 Kick
            particle.position += particle.velocity * dt  # Drift
            particle.force = compute_nbody_force(1, particle, particle.mass,0.5) # Update Acceleration
            particle.velocity = particle.force * dt/2
        t += dt







def lg_fix():
    lg_mass_vals = []
    lg_particle_count = []
    lg_compute_analytic_force = []
    for i in densities[1]:
        if i != 0:
            lg_mass_vals.append(np.log(i))
        else:
            continue

    for i in densities[2]:
        if i != 0:
            lg_particle_count.append(np.log(i))
        else:
            continue

    return lg_mass_vals, lg_particle_count, lg_compute_analytic_force

def draw_figs(number):  # Plots

    if number == 1:
        fig, (ax1,ax2,ax3) = plt.subplots(1,3, figsize =(15,6))
        ax1.plot(lg_fix()[0])
        # plt.yscale('symlog')
        ax1.set_title('Mass Density Function')
        ax1.set_xlabel('Radius')
        ax1.set_ylabel('Mass Density')

        ax2.plot(np.log(densities[0]))
        ax2.set_title('Analytical Hernquist Density Function')
        ax2.set_xlabel('Radius')
        ax2.set_ylabel('Hernquist Density')

        ax3.plot(lg_fix()[1])
        ax3.set_title('Particle Count')
        ax3.set_xlabel('Radius')
        ax3.set_ylabel('Particle Count')
        plt.savefig('HernquistvsMassDensityFunction.png')
        plt.show()

    if number == 2:
        fig, (ax1,ax2,ax3) = plt.subplots(1,3, figsize = (12,8))
        # ax1.plot(sorted(particle_radius),np.log(force_100), label = 'softening = Mass(100)', color = 'blue')
        # ax1.scatter(particle_radius, np.log(force_mis_new), label='softening = Mass(rhm)', color='orange')
        ax1.scatter(particle_radius[0:5000], np.log(force_mis_new*(-1)), label='softening = Mass(rhm)', color='orange')
        ax2.scatter(particle_radius[0:5000], np.log(force_mis_new2*(-1)), label='softening = Mass(rhm)', color='orange')
        ax3.scatter(particle_radius[0:5000], np.log(force_mis_new3*(-1)), label='softening = Mass(rhm)', color='orange')



        # ax2.plot(particle_radius, force_new, label='softening = Mass(rhm)', color='orange')
        # ax3.plot(sorted(particle_radius),force_half, label = 'softening = Mass(rhm/2)', color = 'green')
        ax2.set_xlabel('Particle Index')
        for i in ax1,ax3:
            i.set_xlabel('Radius')
            i.set_ylabel('Force')
        ax1.legend()
        ax2.legend()
        ax3.legend()
        plt.savefig('ForcevsSoftening.png')
        plt.show()

    if number == 3:
        plt.plot(np.log(compute_analytical_force()))
        plt.xlabel('Radius')
        plt.ylabel('Force')
        plt.title('Analytical Force')
        plt.savefig('AnalyticalForce.png')
        plt.show()

    if number == 4:
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
        ax.scatter(particle_position_x, particle_position_y, particle_position_z)
        # plt.scatter(particle_position_x, particle_position_y)
        plt.xlabel('x')
        plt.ylabel('y')
        plt.title('Particle Distribution')
        plt.show()





# compute_nbody_forces(1, softening_mean_inter_sep)
draw_figs(4)
# compute_relaxation()
# leap_frog(2,0.5)



