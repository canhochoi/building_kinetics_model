# -*- coding: utf-8 -*-
"""
Created on Sat May 20 14:08:01 2023

@author: lqluan
"""

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import sympy as sp
import re #for finding regular expression


class ode_construct:    
    def __init__(self, reactant_list, Reaction_list, reactor_eqn, arg):
        self.reactant_list = reactant_list        
        self.reaction_on_catalyst = Reaction_list
        self.reaction_with_reactor = reactor_eqn
        self.parameters = arg
    
    def species_create(reactant_list):
        species = ['S'+str(i) for i in range(len(reactant_list))]
        return species               
   
    def reaction_dict(reaction_list, reactant_list, species):
        """  
    
        Parameters
        ----------
        reaction_list : the list of chemical reactions at each step
        
        species: the list of covnerted chemical species (S0, S1, ...)
        
        reactant_list: the list of chemical species (HCOO*, etc ...)
    
        Returns
        -------
        reactions : a dictionary with keys are each reaction step containing reactants, stoichiometry, products, etc ...
    
        """
        
        reactions = {}
        reaction_convert = ['']*len(reaction_list)
        for i, reaction in enumerate(reaction_list):
            reactants, products = reaction.split(' --> ')
            reactants = reactants.split(' + ')
            products = products.split(' + ')
            
                
            #Step 2a: establish stoichiometry matrix
            reactant_stoichiometries = []
            reactant_orig = []
            for reactant in reactants:        
                    
                #use regular expression package
                stoichiometry, adsorbate = re.match(r'^(\d*)(.*)$', reactant).groups()
                stoichiometry = int(stoichiometry) if stoichiometry else 1
                # if reactant == 'HCOOH(g)':
                #     stoichiometry = Vcat*site_density
                                    
                reactant_stoichiometries.append(stoichiometry)
                reactant_orig.append(adsorbate)
                
            product_stoichiometries = []
            product_orig = []
            for product in products:
                stoichiometry, adsorbate = re.match(r'^(\d*)(.*)$', product).groups()
                stoichiometry = int(stoichiometry) if stoichiometry else 1
                # if product == 'H2(g)' or product == 'CO2(g)':
                #     stoichiometry = Vcat*site_density
                    
                product_stoichiometries.append(stoichiometry)
                product_orig.append(adsorbate)
                
            #Step 2b: convert reactant_orig and product_orig into S species 
            converted_reactants = []
            converted_products = []
            for r in reactant_orig:
                converted_reactants.append(species[reactant_list.index(r)])
            for p in product_orig:
                converted_products.append(species[reactant_list.index(p)])        
                       
            reactions[f"R{i}"] = {
                "reactants": converted_reactants,
                "reactants_stoi": reactant_stoichiometries,
                "products": converted_products,
                "products_stoi": product_stoichiometries,
                "kf_idx": i,
                "kb_idx": i,
                "reactants_orig": reactant_orig,
                "products_orig": product_orig} 
      
        return reactions

    def reactor_dict(reactor_eqn, reactant_list, species, arg_J):
        """
        
    
        Parameters
        ----------
        reactor_eqn : the list of chemical reactions at each step
        
        species: the list of covnerted chemical species (S0, S1, ...)
        
        reactant_list: the list of chemical species (HCOO*, etc ...)
        
        arg_J: arguments of parameters such as reactor volume, R constant, temperature
    
        Returns
        -------
        reactions_J : a dictionary with keys are each reaction step containing reactants, stoichiometry, products, etc ...
        """
        
        V, R, T = arg_J
        reactions_J = {}
        
        for i, reaction_reactor in enumerate(reactor_eqn):
            reactants, products = reaction_reactor.split(' --> ')
            reactants = reactants.split(' + ')
            products = products.split(' + ')    
           
            
            #Step 2a: establish stoichiometry matrix    
            reactant_stoichiometries = []
            reactant_orig = []
            #identify the index of gaseous reactants
            for reactant in reactants:
                if reactant == 'HCOOH_in' or reactant == 'CO2_in' or reactant == 'H2_in':
                    converted_reactants = reactant
                    stoichiometry = np.NaN   #no need to pay attention to
                    reactant_stoichiometries.append(stoichiometry)
                    reactant_orig.append(reactant)
             
            product_stoichiometries = []
            product_orig = []
            for product in products:
                stoichiometry, adsorbate = re.match(r'^(\d*)(.*)$', product).groups()
                # stoichiometry = int(stoichiometry) if stoichiometry else 1
                # stoi is RT/V because the species is gas pressure, not amount of mole
                # dP/dt = RT/V*dn/dt
                # stoichiometry = R*T/V
                stoichiometry = 1 #because species in mole amout of gas
                product_stoichiometries.append(stoichiometry)
                product_orig.append(adsorbate)        
            
            #Step 2b: convert reactant_orig and product_orig into S species 
            converted_reactants = []
            converted_products = []
            for r in reactant_orig:
                if r == 'HCOOH_in' or r == 'CO2_in' or r == 'H2_in':
                    converted_reactants.append(r)
                else:
                    converted_reactants.append(species[reactant_list.index(r)])
            for p in product_orig:
                converted_products.append(species[reactant_list.index(p)])      
            
            reactions_J[f"J{i}"] = {
                "reactants": converted_reactants,
                "reactants_stoi": reactant_stoichiometries,
                "products": converted_products,
                "products_stoi": product_stoichiometries,
                "kf_idx": i,
                "kb_idx": i,
                "reactants_orig": reactant_orig,
                "products_orig": product_orig} 
        return reactions_J
    
    def stoi_matrix(species, n_species, n_reactions, reactions, reactant_list, arg):
        """
        
        Parameters
        ----------
        species: the list of covnerted chemical species (S0, S1, ...)
        n_species : number of chemical species
        n_reactions : number of chemical reactions
        reactions : a dictionary storing all information of chemical reactions
        reactant_list: a list of reactants 
        arg: argument containing parameters such as V, T, R, Vcat, site_density
        Returns
        -------
        stoichiometry matrix of whole reaction

        """
        stoichiometry_matrix = np.zeros((n_species, n_reactions))    

        #can convert into a list
        #list(reactions.values())[0]
        V, R, T, Vcat, site_density = arg
       
        for j, reaction_name in enumerate(reactions.keys()):
            if reaction_name[0] == 'R': 
                reaction = reactions[reaction_name]
                reactants = reaction['reactants']
                products = reaction['products']
                reactant_stoichiometries = reaction['reactants_stoi']
                product_stoichiometries = reaction['products_stoi']
                for i, species_name in enumerate(species):
                    if species_name in reactants:
                        if reactant_list[i] == 'HCOOH(g)':
                            #need Vcat*site_density is needed for amount of mole of HCOOH(g) consumed by catalyst
                            #need RT/V because the species is P_HCOOH = n_HCOOH * RT/V
                            # stoichiometry_matrix[i, j] -= R*T/V*Vcat*site_density 
                            stoichiometry_matrix[i, j] -= Vcat*site_density #because species is the mole amount of gas
                        else:
                            stoichiometry_matrix[i, j] -= reactant_stoichiometries[reactants.index(species_name)]
                    if species_name in products:
                        if reactant_list[i] == 'H2(g)' or reactant_list[i] == 'CO2(g)':
                            # check comments on HCOOH(g)
                            # stoichiometry_matrix[i, j] += R*T/V*Vcat*site_density 
                            stoichiometry_matrix[i, j] += Vcat*site_density       #because species is the mole amount of gas          

                        else:
                            stoichiometry_matrix[i, j] += product_stoichiometries[products.index(species_name)]
            else:
                reaction = reactions[reaction_name]           
                products = reaction['products']
                product_stoichiometries = reaction['products_stoi']
                
                for i, species_name in enumerate(species):
                    if species_name in products:                    
                        stoichiometry_matrix[i, j] += product_stoichiometries[products.index(species_name)]
        
        #apply site conservation
        row_star = reactant_list.index('*')
        stoichiometry_matrix[row_star,:] = -sum(stoichiometry_matrix[0:row_star,:])
        return stoichiometry_matrix

    def get_funcs(n_species):
        """   
    
        Parameters
        ----------
        n_species : number of species reacting 
    
        Returns
        -------
        funcs : time-dependent symbols representing species
    
        """
        t = sp.symbols('t')
        S = sp.symbols('S0:{}'.format(n_species))    
        
        funcs = []
        for i in range(n_species):
            func = sp.Function(f'S{i}')(t)
            funcs.append(func)
        return S, t, funcs
    
    def rate_construct(self, S, funcs, products, product_stoi, klist, k_index):
        """    

        Parameters
        ----------
        S: a tuple of converted chemical species (S0, S1, ...)
        funcs : time dependent species 
        products : chemical symbols of products
        product_stoi : stoichiometry of products
        klist : list of rate constant
        k_index : index of each rate constant 

        Returns
        -------
        rate : net rate of each chemical reaction

        """
        rate = []
        for i, ads in enumerate(products):         
            # Find the index of products[i] in funcs
            S_idx = S.index(sp.Symbol(products[i]))
            if product_stoi[i]: #different from 1
                rate = np.append(rate,funcs[S_idx]**product_stoi[i])
            else:
                rate = np.append(rate,funcs[S_idx])
                
        rate = klist[k_index] * np.prod(rate)    
        return rate

    def reaction_to_equation(self, S, reaction, funcs, kf, kb, arg_J):
        """
        
        Parameters
        ----------
        S: a tuple of converted chemical species (S0, S1, ...)
        reaction : a dictionary storing all information of a reaction 
        funcs : time dependent species 
        kf : forward rate constants
        kb : backward rate constants
        arg_J: arguments of parameters such as reactor volume, R constant, temperature

        Returns
        -------
        net rate of a chemical reaction 

        """
        V, R, T =  arg_J
        reactants = reaction['reactants']
        products = reaction['products']
        kf_idx = reaction['kf_idx']
        kb_idx = reaction['kb_idx']
        reactant_stoi = reaction['reactants_stoi']
        product_stoi = reaction['products_stoi']
        # forward_rate = kf[kf_idx] * np.prod([func for i, func in enumerate(funcs) if reactants[i]])
        # reverse_rate = kb[kb_idx] * np.prod([func for i, func in enumerate(funcs) if products[i]])
        forward_rate = self.rate_construct(S, funcs, reactants, reactant_stoi, kf, kf_idx)
        backward_rate = self.rate_construct(S, funcs, products, product_stoi, kb, kb_idx)
        if 'S7' in reactants:  #has HCOOH(g)
            forward_rate = forward_rate*R*T/V
        if 'S8' in products or 'S9' in products: #has H2(g) or CO2(g)
            backward_rate = backward_rate*R*T/V
        return forward_rate - backward_rate
    
    def reactor_to_equation(S, reaction, funcs, kf_J, kb_J, n_in, v_in, v_out, arg_J):
        """
        
        Parameters
        ----------
        S: a tuple of converted chemical species (S0, S1, ...)
        reaction : a dictionary storing all information of a reaction 
        funcs : time dependent species 
        kf : forward rate constants
        kb : backward rate constants
        arg_J: arguments of parameters such as reactor volume, R constant, temperature
        
        Returns
        -------
        net rate of a chemical reaction 

        """
        V, R, T = arg_J
        reactants = reaction['reactants']
        products = reaction['products']
        kf_idx = reaction['kf_idx']
        kb_idx = reaction['kb_idx']
        reactant_stoi = reaction['reactants_stoi']
        product_stoi = reaction['products_stoi']
        # forward_rate = kf[kf_idx] * np.prod([func for i, func in enumerate(funcs) if reactants[i]])
        # reverse_rate = kb[kb_idx] * np.prod([func for i, func in enumerate(funcs) if products[i]])
        S_idx = S.index(sp.Symbol(products[0]))    
          
        # multiply RT/V because the species is gas pressure, not amount of mole
        # r = n_in[reactants[0]]*v_in[reactants[0]] - funcs[S_idx]*V/(R*T)*v_out[reactants[0]]/V
        r = n_in[reactants[0]]*v_in[reactants[0]] - funcs[S_idx]*v_out[reactants[0]]/V

        return r
    
    def ode_build(species, ode_matrix, t, funcs):
        """
        
        Parameters
        ----------
        species : list of chemical species
        ode_matrix : stoichiometry matrix * reaction list 
        t: a symbolic converted symbol for time 
        funcs: an array of time dependent species 
        Returns
        -------
        equations : time-dependent equations for ode solver
        """
        eqs = []

        for i, equation in enumerate(ode_matrix):
            eq = sp.Eq(sp.Function(species[i])(t).diff(t), equation)
            eqs.append(eq)

        # Convert the equations to Python functions
        equations = [sp.lambdify((t, *funcs), eq.rhs) for eq in eqs]

        return equations