# -*- coding: utf-8 -*-
"""
Created on Mon Jul  4 17:53:26 2022
@author: op5
Code for paper "Multifarious selection on the evolution of visual signals",
<add ref of paper>

The code computes the temporal evolution of the distribution of phenotypes given 
random mutations and three factors that shape the fitness landscape, namely the 
cost of each phenotype, the forces of sexual selection and the forces of natural
selection. The evolution of the whole population's distribution of phenotypes is 
described using a Fokker-Planck equation. In this equation, the mutations are 
modelled using the diffusion term and the fitness landscape is modelled by the 
free energy term (accordingly, maximizing the population fitness is equivalent
to minimizing the population's free energy, ie "fitness = -free energy"; see text
for relationship between free energy and Fokker-Planck equation)  

NB1: The code mainly uses functions created by Amanuel Wolde-Kidan,
https://github.com/woldeaman/numerical_2D_FPE_solver (retrieved June 2022)
Please cite Wolde-Kidan, A. et al. Particle Diffusivity and Free-Energy Profiles
in Hydrogels from Time-Resolved Penetration Data. Biophys. J. 120, 463-475, 
doi:10.1016/j.bpj.2020.12.020 (2021) and <add ref of paper> if you use this code.

NB2: The free energy components that define the total energy landscape F= 100*Ztot
are defined in an annex Matlab function called create_landscape_for_scenarios_ABC.m

"""

# -*- coding: utf-8 -*-
# first start for fokker-planck equation in 2D
import numpy as np
import time
import sys
import matplotlib.pyplot as plt
import scipy.linalg as al
import numpy.linalg as la
from mpl_toolkits.mplot3d.axes3d import Axes3D  # needed for 3D plot support
from matplotlib import cm
from mat4py import loadmat # needed to get import Matlab data
import os
import scipy.io # needed to export output in Matlab *.mat format (requires numpy as well)

# This will be the output directory
if not os.path.exists('output_test'):
    os.makedirs('output_test')

startTime = time.time()

#=================================     
#  choose scenario to implement
#=================================   
# please comment/uncomment manually
#scenario = 'scenarioA';
#scenario = 'scenarioB1';
#scenario = 'scenarioB2'; # NB: not really interesting; we just swap the weak (Fsex) and strong forces (Fnat) wrt 'scenarioB1'
scenario = 'scenarioC';


def plotting(X, Y, F, D, tt, cc, savePath=''):
    # plotting results
    fig = plt.figure()  # D and F in one figure
    ax = fig.add_subplot(121, projection='3d')
    ax.plot_surface(X, Y, F, cmap=cm.coolwarm, antialiased=True)
    # start adding
#    dx, dy = np.gradient(F)
#    ax.quiver(X,Y,X,Y,scale=5,angles="uv",headwidth=5)
#    ax.quiver(X,Y,dx,dy,scale=5,angles="uv",headwidth=5)
#    ax.quiver(X, Y, dx, dy)#,scale=5,angles="uv",headwidth = 5)       
    # end adding        
    ax.set_title('Free Energy')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('F [k$_B$T]')
    ax.view_init(90, -90) # added op; format is (elevation, azimuth)
    ax = fig.add_subplot(122, projection='3d')
    ax.plot_surface(X, Y, D, cmap=cm.coolwarm, antialiased=True)
    ax.set_title('Diffusivity')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('D [bins$^2$/timestep]')
    ax.view_init(90, -90) # added op; format is (elevation, azimuth)
    fig.suptitle( 'Landscapes for ' + scenario)
    plt.savefig(savePath+'output_v2\\'+'DF_' + scenario + '_N50_start15_15.png', dpi = 600)
    #plt.savefig(savePath+'output_v2\\'+'DF_' + scenario + '_N50_start20_12.png', dpi = 600)

    cMax = np.max(cc)
    cMin = np.min(cc)

    M = tt.size
    fig = plt.figure()  # concentration profiles in second figures
    fig.set_size_inches(12, 11) # added op
    row = np.ceil(M/2)  # number of rows for c-plot
    col = np.floor(M/2)  # number of columns for c-plot

    for i in range(M):
        ax = fig.add_subplot(row, col, i+1, projection='3d')
        # ax.plot_surface(X, Y, cc[i], cmap=cm.coolwarm, vmin=cMin, vmax=cMax/10, label='Original', antialiased=True) # original code
        ax.plot_surface(X, Y, cc[i], cmap=cm.coolwarm, vmin=cMin, vmax=np.max(cc[i]), label='Original', antialiased=True)
        ax.set_title('t = %i' % tt[i])
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Count')
        ax.set_zlim(cMin, np.max(cc[i]))  # was ax.set_zlim(cMin, cMax)
        ax.view_init(90, -90) # added op; format is (elevation, azimuth)
    fig.suptitle( 'Population evolution for ' + scenario)    
    plt.savefig(savePath+'output_v2\\'+'profiles_' + scenario + '_N50_start15_15.png', dpi = 600)
    #plt.savefig(savePath+'output_v2\\'+'profiles_' + scenario + '_N50_start20_12.png', dpi = 600)


def WMatrix(D, F, dimY=None, dx=1):
    '''
    Function computes rate matrix for 2D systems, here the reduced index i
    is used defined as i = dimY*x + y, with dimY - number of bins in y-disc.
    D and F have to be reduced to 1D-index too in order for this to work!
    '''
    # using auxillary functions gFunc and fFunc as proposed by
    # Grima et. al - "Accurate discretization of advection-diffusion equations"
    # (2004) PRE https://doi.org/10.1103/PhysRevE.70.036703
    def g(i):
        return np.sqrt(D[i])*np.exp(F[i]/2)  # F in kB*T

    def f(i):
        return np.sqrt(D[i])*np.exp(-F[i]/2)  # F in kB*T
    # wRate(i --> j) = (f(j)*g(i))/(dx**2)

    N = D.size  # size of WMatrix is (dimX*dimY)x(dimX*dimY)
    if dimY is None:
        dimY = np.sqrt(N)  # equal discretization assumed if not specified

    #  first compute transition rates to come to bin i
    W = np.array([[f(i-1)*g(j-1)/(dx**2)
                   # only possible to go to neighbouring bins
                   if(
                       # nearest neighbours for first column
                       ((i-1) % dimY == 0 and
                        (abs(j-i) == dimY or
                         j == i+1 or
                         j == i+dimY+1 or
                         j == i-dimY+1))
                       or
                       # nearest neighbours for last column
                       (i % dimY == 0 and
                        (abs(j-i) == dimY or
                         j == i-1 or
                         j == i+dimY-1 or
                         j == i-dimY-1))
                       or
                       # nearest neighbours central bins
                       (((i-1) % dimY != 0 and i % dimY != 0) and
                        (abs(j-i) == dimY or
                         abs(j-i) == 1 or
                         abs(j-i-dimY) == 1 or
                         abs(j-i+dimY) == 1))
                       )
                   else 0
                   for j in range(1, N+1)] for i in range(1, N+1)])
    # indices are shifted because of modulo operator, doesn't handle 0 properly

    # then add rates to leave on main diagonal from original matrix
    for i in range(N):
        W[i, i] = -np.sum([W[j, i] for j in range(N)])

    if np.any(np.sum(W, 0) > 1E-10):
        print('WMatrix not row stochastic!')

    return W


def computeC(cc, tt, W=None, T=None, D=None, F=None, dimY=None, dx=1):
    '''
    Calculates concentration profiles at time t from W or T matrix with
    reflective boundaries, based on concentration profile cc
    '''
    # calculate only variables that are not given
    if T is None:
        if W is None:
            if (D is None) or (F is None):
                print('Error: You gotta give me something... no T, W, D or '
                      'F given!')
                sys.exit()
            else:
                W = WMatrix(D, F, dimY=dimY, dx=dx)
                T = al.expm(W)  # exponential of W
        else:
            T = al.expm(W)  # exponential of W

    return np.dot(la.matrix_power(T, tt), cc)


def main():
    # |- - - - - - > y             layout is in a way that y-axis points to
    # |(0, 0) (0, 1) ...           the right of matrix and x-axis points down,
    # |(1, 0) ...                  origin (0, 0) is in upper left corner,
    # | ...                        like a coordinate system tilted
    # v x                          at a right angle

    # defining grid and initial distribution (should match the size of Ztot, imported from Matlab!)
    dimX = 50
    dimY = 50
    cInit = 100
    c0 = np.zeros((dimX, dimY))
    # define initial population
    startX = 15 # default is 15/ alternative startX for scenarioC is 20
    startY = 15 # default is 15/ alternative startY for scenarioC is 12
    c0[startX-1:startX+1, startY-1:startY+1] = cInit
    # meshgrid function needs X, Y in switchted order for correct ouput
    X, Y = np.meshgrid(np.arange(dimY), np.arange(dimX))

    # =====================================
    #   setting D and F
    # =====================================
    # diffusion (does not depend on scenario)
    D = np.ones((dimX, dimY))/(0.75*100)  # flat diffusivity profile (was 100 originally)
    
    # all forces in free energy term (scenarios created with Matlab, see material_polymorphism_paper_v3.m)
    if(scenario == 'scenarioA'):
        data = loadmat('total_free_energy_landscape_for_Fokker_Plank_equation_scenarioA_Siz50.mat')
    elif(scenario == 'scenarioB1'):
        data = loadmat('total_free_energy_landscape_for_Fokker_Plank_equation_scenarioB1_Siz50.mat')
    elif(scenario == 'scenarioB2'):
        data = loadmat('total_free_energy_landscape_for_Fokker_Plank_equation_scenarioB2_Siz50.mat')
    elif(scenario == 'scenarioC'):
        data = loadmat('total_free_energy_landscape_for_Fokker_Plank_equation_scenarioC_Siz50.mat')
    F = np.array( data["Ztot"] )
    F = 100*F
    #======================================     
    #           end of scenarios
    #======================================       

    # ===== compute profiles from given diffusion D and free energy landscape =====
    # a few time points only
    # tt = np.array([0, 100, 200, 500, 750, 1000])
    # standard implementation (ie, "medium" time resolution)
    tt = np.linspace(0, 1000, 41, dtype=np.integer)
    # higher temporal resolution (e.g., for creating videos of the dynamic of the distribution of the population phenotypes, see Suppl. Material)
    # tt = np.linspace(0, 1000, 81, dtype=np.integer)

    print(tt)
    cInput = [computeC(c0.reshape(c0.size), tt[i], D=D.reshape(D.size), F=F.reshape(F.size),
                       dimY=dimY).reshape(dimX, dimY) for i in range(tt.size)]
    print([np.sum(cInput[i]) for i in range(len(cInput))])  # check for conservation of mass
    
    # show the starting distribution of phenotypes
    fig = plt.figure()  # D and F in one figure
    ax = fig.add_subplot(121, projection='3d')
    ax.plot_surface(X, Y, cInput[0], cmap=cm.coolwarm, antialiased=True)
    
    
    # ===== standard implementation =====
    if(scenario == 'scenarioA') or (scenario == 'scenarioB1') or (scenario == 'scenarioB2'):
        name_out = 'cInput_across_time_ie_population_distribution_evolution_' + scenario + '.mat'
    elif(scenario == 'scenarioC'): # comment/uncomment manually
        name_out = 'cInput_across_time_ie_population_distribution_evolution_' + scenario + '_15_15' + '.mat'
        # name_out = 'cInput_across_time_ie_population_distribution_evolution_' + scenario + '_12_20' + '.mat'
        
    # ===== for higher temporal resolution =====
#    if(scenario == 'scenarioA') or (scenario == 'scenarioB1') or (scenario == 'scenarioB2'):
#        name_out = 'cInput_across_time_ie_population_distribution_evolution_' + scenario + '_highTempRes.mat'
#    elif(scenario == 'scenarioC'): # comment/uncomment manually
#        name_out = 'cInput_across_time_ie_population_distribution_evolution_' + scenario + '_15_15' + '_highTempRes.mat'
#        #name_out = 'cInput_across_time_ie_population_distribution_evolution_' + scenario + '_12_20' + '_highTempRes.mat'
       
    # export cInput in Matlab format
    scipy.io.savemat(name_out, ({'cInput' : cInput}))
    
    # start adding (uncomment if needed)
    #    dx, dy = np.gradient(F)
    #    ax.quiver(X,Y,X,Y,scale=5,angles="uv",headwidth=5)
    #    ax.quiver(X,Y,dx,dy,scale=5,angles="uv",headwidth=5)
    #    ax.quiver(X, Y, dx, dy)#,scale=5,angles="uv",headwidth = 5)       
    # end adding        
    
    ax.set_title('Start. distr. of phenotypes')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Relative freq.')
    ax.view_init(90, -90) # added op; format is (elevation, azimuth)
    ax = fig.add_subplot(122, projection='3d')
    ax.plot_surface(X, Y, cInput[5], cmap=cm.coolwarm, antialiased=True)
    ax.set_title('Distr. of phenotypes, intermed. step')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Relative freq.')
    ax.view_init(90, -90) # added op; format is (elevation, azimuth)

    # plot output (comment when generating long chains)
    # plotting(X, Y, F, D, tt, cInput)

if __name__ == "__main__":
    main()
    print("Execution time was %.2f minutes"
          % ((time.time() - startTime)/60))

