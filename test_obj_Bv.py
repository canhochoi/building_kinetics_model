# -*- coding: utf-8 -*-
"""
Created on Sun May 21 16:11:22 2023

@author: lqlua
"""
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import sympy as sp
from ode_obj import ode_construct

T = 300 + 273.15
R = 8.314 #J/mol.K
NA = 6.023*10**23 #atoms per mole

n_HCOOH = 1 #mol/L
n_H2 = 0
n_CO2 = 0

vin_HCOOH = 7e-1   #L/s
vout_HCOOH = 7.1e-1
vin_H2 = 0
vout_H2 = vout_HCOOH
vin_CO2 = 0
vout_CO2 = vout_HCOOH

v_in = {'HCOOH_in': vin_HCOOH,
        'H2_in': vin_H2,
        'CO2_in':vin_CO2}
v_out = {'HCOOH_in': vout_HCOOH,
        'H2_in': vout_H2,
        'CO2_in':vout_CO2}

n_in = {'HCOOH_in': n_HCOOH,
        'H2_in': n_H2,
        'CO2_in': n_CO2}

rhosite = 1
Sa = 1000 #m2/g
# Sa = 70e6 #m2/g #catalyst surface
rho_c = 2e6 #g/m3  #catalyst density
Asite = 1e-19 #  m^2/site   #area of a site
phi = 0.6   #porousity 
V = 0.8e-4  #L volume of reactor #conversion_pc = 0.00025
Vcat = 1e-6*210  #to compensate for reduce in surface area 
site_density = Sa*rho_c/Asite*(1-phi)/1000   #site/L 
site_density = site_density/NA    #mole/L


reactant_list = ['HCOOH*', 'HCOObi*', 'HCOOmono*', 'CO2*', 'H*', 'H2*', '*', 'HCOOH(g)', 'CO2(g)', 'H2(g)']

Reaction_list = ['HCOOH(g) + * --> HCOOH*',
                 'HCOOH* + * --> HCOObi* + H*',
                 'HCOObi* --> HCOOmono*',
                 'HCOOmono* + * --> H* + CO2*',
                 '2H* --> H2* + *',
                 'CO2* --> CO2(g) + *',
                 'H2* --> H2(g) + *',
                 'HCOOH* + * --> HCOOmono* + H*']

reactor_eqn = ['HCOOH_in --> HCOOH(g)',               
               'CO2_in --> CO2(g)',              
               'H2_in --> H2(g)']


kf_J = [vin_HCOOH, vin_CO2, vin_H2]
kb_J = [vout_HCOOH, vout_CO2, vout_H2]

kf_list = np.array([1.62241997e+03, 1.19801121e+13, 3.79678473e+07, 3.61031286e+06,
       4.34096330e-01, 6.00154146e+11, 3.14145729e+10, 1.58177601e+13])

kb_list = np.array([4.02889628e+03, 7.77506084e+02, 2.54703711e+10, 3.87451827e+03,
       9.75158147e+12, 1.65923571e+03, 7.75267240e+03, 6.88663707e+05])

n_species = len(reactant_list)
n_reactions = len(Reaction_list) + len(reactor_eqn)

arg = (V, R, T, Vcat, site_density)


#output the object
cat_obj = ode_construct(reactant_list, Reaction_list, reactor_eqn, arg)


# Step 1: Convert reactants to 'S0', 'S1', 'S2', 'S3'
species = ode_construct.species_create(reactant_list)

# Step 2: Build reactions dictionary

#reactions on catalyst
reactions = ode_construct.reaction_dict(Reaction_list, reactant_list, species)

#reactions involved reactor
arg_J = (V, R, T)
reactions_J = ode_construct.reactor_dict(reactor_eqn, reactant_list, species, arg_J)

#update reactions_J to reactions
reactions.update(reactions_J)

# Step 3: Construct reaction stoichiometry matrix for net rate

stoichiometry_matrix = ode_construct.stoi_matrix(species, n_species, n_reactions, reactions, reactant_list, arg)

# Step 4: Define the symbols and functions
S, t, funcs = ode_construct.get_funcs(n_species)

#Step 5: Define array of net rate for each reaction

# rate_reaction = ode_construct.get_net_rate(reactions, S, funcs)

rate_reaction = []
for reaction_name in reactions.keys():    
    if reaction_name[0] == 'R': 
        reaction = reactions[reaction_name]
        r = cat_obj.reaction_to_equation(S, reaction, funcs, kf_list, kb_list, arg_J)
    else:
        reaction = reactions[reaction_name]
        r =  ode_construct.reactor_to_equation(S, reaction, funcs, kf_J, kb_J, n_in, v_in, v_out, arg_J)
        
    rate_reaction.append(r)
    
# convert the expressions to a numpy array
rate_reaction = np.array(rate_reaction)

#Step 6: obtain ode matrix by multiplying the stoichiometry matrix by the expressions array
ode_matrix = np.dot(stoichiometry_matrix, rate_reaction)


#Step 7: define system of differential equations

equations = ode_construct.ode_build(species, ode_matrix, t, funcs)

# Define the initial conditions
t0 = 0
idx = reactant_list.index('*')
y0 = [0]*len(reactant_list)
y0[idx] = rhosite

# Define the time interval
t_span = [0, 8000]


# Solve the system of differential equations
sol = solve_ivp(lambda t, vars: [equation(t, *vars) for equation in equations],
                t_span, y0, method = 'Radau', atol = 1e-25, rtol = 1e-10)

time_sol = sol.t
y_sol = sol.y

#0 - 1 - 2: HCOOH* - HCOObi* - HCOOmono*
#3 - 4 - 5: CO2* - H*  - H2*      
#6 - 7 - 8: * - HCOOH - CO2
#9        : H2 (mole)
dct = {}
for i, reactant in enumerate(reactant_list):
    dct[reactant] = y_sol[i,:] 

lines = []   

#coverage at steady state
coverage = []

fig = plt.figure()
ax = fig.add_subplot(111)
# for i,name_species in enumerate(species):
for name_species in list(dct.keys()):    
    if name_species != 'HCOOH(g)' and name_species != 'H2(g)' and name_species != 'CO2(g)':
        coverage.append(dct[name_species][-1]/rhosite)
        lines += plt.plot(time_sol,dct[name_species]/rhosite, label = name_species)

ax.set_xlabel('Time (s)')
ax.set_ylabel('Coverage')
labels = [l.get_label() for l in lines]
ax.legend(lines, labels)
    
fig = plt.figure()
ax = fig.add_subplot(111)
labels = []
for name_species in list(dct.keys()): 
    if name_species == 'HCOOH(g)' or name_species == 'H2(g)' or name_species == 'CO2(g)':
        plt.plot(time_sol,dct[name_species])
        labels.append(name_species)
ax.set_xlabel('Time (s)')
ax.set_ylabel('Pressure (Pa)')
ax.legend(labels)

fig = plt.figure()
ax = fig.add_subplot(111)
plt.bar(reactant_list[0:-3],coverage)

#check flux balance
f = []
for rate in rate_reaction:
    f.append(sp.lambdify((t, *funcs), rate ))


t_idx = -1
val = []    
for rate in f:
    val.append(rate(t_idx, *y_sol[:,t_idx]))

ratio = np.zeros((6,1))
#r0 = r1 + r7
ratio[0] = val[0]/(val[1]+val[7])
#r1 = r2 
ratio[1] = val[1]/(val[2])
#r3 = r0 
ratio[2] = val[3]/(val[0])   
#r5 = r3 
ratio[3] = val[5]/(val[3])
#2r4 = r1 + r3 + r7
ratio[4] = 2*val[4]/(val[1] + val[3] + val[7])
#r4 = r6
ratio[5] = val[4]/val[6]
