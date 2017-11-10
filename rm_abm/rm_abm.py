from mesa import Model, Agent
from mesa.datacollection import DataCollector
from .wolfsheep_schedule import RandomActivationByBreed # from WolfSheep
from .evolvable import EvolvableVector
from matplotlib import pyplot as plt
from collections import defaultdict
from intervaltree import Interval, IntervalTree
import numpy as np
import pandas as pd
import random

def by_genotype(g):
    '''Return a lambda to filter by genotype'''
    def wrapper(a):
        if a.genotype == g:
            return True
        else:
            return False
    return wrapper

def by_methylation(m):
    '''Return a lambda to filter by methylation'''
    def wrapper(a):
        if a.methylation == m:
            return True
        else:
            return False
    return wrapper

def bernoulli(p):
    '''Return a bernoulli random variable'''
    return np.random.binomial(1,p)

def index_of_first(lst, pred):
    '''Get the index of the first object where pred is true'''
    for i,v in enumerate(lst):
        if pred(v):
            return i
    return None

def get_breed_filtered_count(breed_class, L):
    '''
    Returns the current number of agents of a certain breed
    where L -> True
    '''
    def wrapper(model):
        return len(list(filter(L, model.schedule.agents_by_breed[breed_class].values())))
    return wrapper


class BaseModel(Model):
    '''
    Phage-Bacteria with RM systems
    '''

    def __init__(self,
                 initial_phage = 10,
                 initial_fraction_p_m1 = 1,
                 initial_fraction_p_g1 =1,
                 initial_bacteria=100,
                 fraction_b_m1=0.5,                 
                 phage_inactivation_time=1,
                 phage_burst_size=3,
                 phage_off_diagonal=0.01,
                 phage_mutation_step = 0,
                 phage_mutation_freq = 0,
                 re_degrade_foreign_0 = 1e-3,
                 re_degrade_foreign_1 = 1e-3,
                 bacteria_per_step = 10,
                 encounter_width = 0.01,
                 verbose=False,
                 latency = 0.5,
                 epi_inheritance = 1):
        '''
        Create a new Phage-Bacteria model with the given parameters.

        Args:
        initial_phage (int)               number of phage to start with
        intial_fraction_p_rm1 (float)     percent of phage with methylation 1
        initial_fration_p_g1 (float)      percent of phage with genotype 1
        initial_bacteria (int)            number of bacteria to start with
        fraction_b_m1 (float)             the fraction of bacteria with R-M system 1 that 
                                          are added each epoch
        phage_inactivation_time (int)     number of epochs a phage can live outside a host
        phage_burst_size (int)            number of phage produced when it lyses a cell
        phage_off_diagonal (float)        starting affinity of (e.g.) phage genotype 1 for
                                          bacteria genotype 1
        phage_mutation_step (float)       how big mutations are when they happen
        phage_mutation_freq (float)       probability of a mutation 
        re_degrade_foreign (float)        probability that a R-M system will not kill a phage 
                                          with the opposite methylation pattern
        bacteria_per_step (int)           number of bacteria added during each epoch
        encounter_width (int)             the world is [0,1], encounter width is the amount of 
                                          bacteria that can be encountered by the phage
        verbose (bool)                    print info about the simulation
        latency (float)             
        epi_inheritance (float or string) odds that the progeny gets the parent's methylation state. Otherwise gets None. If "genetic" then the phage inherits its parent's methylation state
        '''
        
        # set parameters
        self.initial_phage = initial_phage
        self.initial_bacteria = initial_bacteria
        self.initial_fraction_p_m1 = initial_fraction_p_m1
        self.initial_fraction_p_g1 = initial_fraction_p_g1
        self.fraction_b_m1 = fraction_b_m1
        self.verbose = verbose
        self.bacteria_per_step = bacteria_per_step
        self.phage_off_diagonal = phage_off_diagonal
        self.re_degrade_foreign_0 = re_degrade_foreign_0
        self.re_degrade_foreign_1 = re_degrade_foreign_1
        self.phage_mutation_step = phage_mutation_step
        self.phage_mutation_freq = phage_mutation_freq
        self.phage_inactivation_time = phage_inactivation_time
        self.phage_burst_size = phage_burst_size
        self.encounter_width = encounter_width
        self.agent_width = 0.0001
        self.latency = latency
        self.epi_inheritance = epi_inheritance

        if self.encounter_width > 1 or self.encounter_width < 0:
            raise ValueError("Encounter width must be between 0 and 1")
        
        self.schedule = RandomActivationByBreed(self)
        
        model_reporters = {
            "phage" : lambda m : m.schedule.get_breed_count(Phage),
            "bacteria" : lambda m :  m.schedule.get_breed_count(Bacteria),
            "bacteria_meth_0" : lambda m: get_breed_filtered_count(Bacteria,by_methylation(0))(m),
            "phage_meth_0" : lambda m: get_breed_filtered_count(Phage,by_methylation(0))(m)
        }

        agent_reporters = {
            "breed" : lambda a : a.breed,
            "methylation" : lambda a: a.methylation,
            "genotype" : lambda a: a.genotype,
            "inactivation" : lambda a: a.inactivation}
        
        self.datacollector = DataCollector(model_reporters = model_reporters,
                                           agent_reporters=agent_reporters)

        self.current_ID = 0
        # Create phage
        self.add_phage()
        #Create bacteria
        self.add_bacteria(self.initial_bacteria)

        self.running = True

    def get_next_ID(self):
        self.current_ID += 1
        return self.current_ID

    def add_phage(self):
        #Create phage
        #phage start with affinity according to a symmetric matrix
        affinity = np.array([[1-self.phage_off_diagonal, self.phage_off_diagonal],
                             [self.phage_off_diagonal, 1-self.phage_off_diagonal]])
        
        for i in range(self.initial_phage):
            # sample methylation
            rm_probs = [1-self.initial_fraction_p_m1,
                        self.initial_fraction_p_m1]
            rm = np.random.choice([0,1],p=rm_probs)
            # sample genotypes
            g_probs = [1-self.initial_fraction_p_g1,
                       self.initial_fraction_p_g1]
            g = np.random.choice([0,1], p=g_probs)
            # assign affinity based on genotype
            p_affinity = EvolvableVector(
                affinity[g,:].copy(),
                self.phage_mutation_step,
                self.phage_mutation_freq)
            phage = Phage(
                self,
                self.get_next_ID(),
                g,
                rm,
                self.phage_inactivation_time,
                p_affinity,
                0) # parent is 0 for first generation
            self.schedule.add(phage)

    def add_bacteria(self, num):
        for i in range(num):
            g = np.random.choice([0,1],p=[1-self.fraction_b_m1, self.fraction_b_m1])
            if g == 0:
                bacteria = Bacteria(self, self.get_next_ID(), g, g,
                                    self.re_degrade_foreign_0)
            elif g == 1:
                bacteria = Bacteria(self, self.get_next_ID(), g, g,
                                    self.re_degrade_foreign_1)
            else:
                raise(ValueError("Unknown genotype"))
            self.schedule.add(bacteria)
            
    def step(self):
        self.datacollector.collect(self)
        # Shuffle the location of agents each time
        self.tree = IntervalTree()
        for unique_id, agent in self.schedule.agents_by_breed[Bacteria].items():
            pos = random.random()
            self.tree.add(Interval(pos, pos + self.agent_width, agent))
        self.schedule.step()
        if self.verbose:
            print([self.schedule.time,
                   self.schedule.get_breed_count(Phage),
                   self.schedule.get_breed_count(Bacteria)])
        self.add_bacteria(self.bacteria_per_step)
        
    def run_model(self, step_count=200):
        if self.verbose:
            print('Initial number phage: ', 
                self.schedule.get_breed_count(Phage))
            print('Initial number bacteria: ', 
                self.schedule.get_breed_count(Bacteria))

        for i in range(step_count):
            self.step()

        if self.verbose:
            print('')
            print('Final number phage: ', 
                self.schedule.get_breed_count(Phage))
            print('Final number bacteria: ',
                self.schedule.get_breed_count(Bacteria))


class SpikeIn(BaseModel):

    def __init__(self,
                 spike_in_affinity_0 = 0,
                 spike_in_methylation = 0,
                 **kwargs):

        BaseModel.__init__(self, **kwargs)

        self.spike_in_affinity_0 = spike_in_affinity_0
        self.spike_in_methylation = spike_in_methylation

        p_affinity = EvolvableVector(np.array([spike_in_affinity_0, 1-spike_in_affinity_0]),
                                     self.phage_mutation_step,
                                     self.phage_mutation_freq)

        for i in range(1,self.phage_burst_size,1):
            phage = Phage(
                self,
                i*-10, 
                0, #genotype 0
                self.spike_in_methylation,
                self.phage_inactivation_time,
                p_affinity,
                -1) #all have parent -1
        
            self.schedule.add(phage)
        
        


class Phage(Agent):
    '''A phage that recognizes and infects bacteria'''
    def __init__(self, model, unique_id, genotype,
                 methylation, inactivation,
                 affinity, parent, last_infected = None):
        '''
        Args:
        affinity : an EvolvableVector with indicies corresponding 
                   to bacterial genotypes
        '''
        self.model = model
        self.genotype = genotype
        self.methylation = methylation
        self.inactivation = inactivation
        self.affinity = affinity
        self.unique_id = unique_id
        self.breed = "Phage"
        self.parent = parent
        self.dead = False
        self.last_infected = last_infected
        if not unique_id:
            print(unique_id)
            raise ValueError

    def step(self):
        inactivated = self.inactivate()
        if not inactivated:
            self.infect()
        
    def inactivate(self):
        '''Remove from schedule if longer than inactivation parameter. Or if
        phage infected last step.  Else decrement counter.

        '''
        self.inactivation -= 1
        if (self.inactivation < 0) or self.dead:
            self.model.schedule.remove(self)
            return True
        
    def infect(self):
        ''' '''
        pos = random.random()
        encounter_radius = self.model.encounter_width/2
        if (pos - encounter_radius) < 0: # if near lower boundary
            remainder = abs(pos - encounter_radius)
            lwr = self.model.tree[0:(pos + encounter_radius)]
            upr = self.model.tree[(1 - remainder):1]
            encountered = [i.data for i in lwr.union(upr)]
        elif (pos + encounter_radius) > 1:
            remainder = (pos + encounter_radius) - 1
            lwr = self.model.tree[0:remainder]
            upr = self.model.tree[(pos - encounter_radius):1]
            encountered = [i.data for i in lwr.union(upr)]
        else:
            intervals = self.model.tree[(pos - encounter_radius):\
                                  (pos + encounter_radius)]
            encountered = [i.data for i in intervals]
        if len(encountered) > 0:
            random.shuffle(encountered)
            genotypes = map(lambda x: x.genotype, encountered)
            genotypes = np.array(list(genotypes))
            infection_probs = self.affinity.vector[genotypes]
            ## Infection probability = odds of infecting a genotype if
            ## it's the only bacteria in the encounter
            ## radius. Otherwise its infection probability/n where n
            ## is the number of agents
            #if np.sum(infection_probs) > 0:
            bernoullis = map(bernoulli, infection_probs)
            to_infect = index_of_first(bernoullis, lambda x: x == 1)
            if to_infect is not None:
                target = encountered[to_infect]
                if target.phage is not None: # die if already infected
                    self.dead = True
                else:
                    target.phage = self
                    self.dead = True # set to remove from schedule

class Bacteria(Agent):
    '''A bacteria with a methylation pattern and a coat protein type'''
    
    def __init__(self, model, unique_id, genotype, methylation, re_degrade_foreign):
        self.model = model
        self.genotype = genotype
        self.methylation = methylation
        self.phage = None
        self.re_degrade_foreign = re_degrade_foreign
        self.unique_id = unique_id
        self.breed = "Bacteria"
        self.inactivation = np.nan
        self.last_infected = None # This just simplifies the agent reporters.
        self.established = False
        if not unique_id:
            raise ValueError
    
    def step(self):
        if self.phage:
            if self.maybe_degrade():
                ## phage dies
                self.phage = None
                self.infected = None
            elif bernoulli(self.model.latency):
                ## cell dies, phage are produced
                if (self.model.epi_inheritance != -1) and\
                   (self.model.epi_inheritance != -2): # This is epigenetic inheritance
                    for i in range(self.model.phage_burst_size):
                        if bernoulli(self.model.epi_inheritance):
                            methylation = self.methylation # Inherits from bacteria
                        else:
                            methylation = None # Can get a "None" inheritance
                        phage = Phage(
                            self.model,
                            self.model.get_next_ID(),
                            self.phage.genotype,
                            methylation,
                            self.model.phage_inactivation_time,
                            self.phage.affinity.copy(), ## mutations possible here
                            self.phage.unique_id,
                            last_infected = self.genotype) 
                        self.model.schedule.add(phage)
                elif self.model.epi_inheritance == -1: # This is genetic inheritance
                    for i in range(self.model.phage_burst_size):
                        phage = Phage(
                            self.model,
                            self.model.get_next_ID(),
                            self.phage.genotype,
                            self.phage.methylation, #inherits parent's methylation,
                            self.model.phage_inactivation_time,
                            self.phage.affinity.copy(), ## mutations possible here
                            self.phage.unique_id,
                            last_infected = self.genotype) 
                        self.model.schedule.add(phage)
                elif self.model.epi_inheritance == -2: #Flip a coin for methylations state
                    for i in range(self.model.phage_burst_size):
                        phage = Phage(
                            self.model,
                            self.model.get_next_ID(),
                            self.phage.genotype,
                            np.random.choice([0,1]), #Random
                            self.model.phage_inactivation_time,
                            self.phage.affinity.copy(), ## mutations possible here
                            self.phage.unique_id,
                            last_infected = self.genotype) 
                        self.model.schedule.add(phage)
                else:
                    raise ValueError("Don't know how to deal with inheritance")    
                self.model.schedule.remove(self)


    def maybe_degrade(self):
        '''Degrades phage according to a bernoulli flip'''
        if not self.established:
            if self.methylation == self.phage.methylation:
                return False
            else:
                if bernoulli(self.re_degrade_foreign):
                    return True
                else:
                    self.established = True
                    return False
        else:
            return False
