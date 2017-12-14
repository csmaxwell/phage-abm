import numpy as np

class EvolvableVector(object):
    '''A class to implement an evolving vector. Returns a mutated
    version of the vector when copied'''

    def __init__(self, vector,  mutation_size, mutation_frequency):
        '''
        Args:
        
        vector (np.array): the vector that should be evolving
        mutation_size (float): size of mutational jumps
        mutation_frequency (float): probability of a jump
        '''
        
        self.vector = self.constrain(vector)
        self.mutation_size = mutation_size
        self.mutation_frequency = mutation_frequency

    def __str__(self):
        return self.vector.__str__()
    
    def mutate(self, vector):
        '''
        Step through the vector, make mutations to each entry.
        '''
        def get_sign(x):
            if x:
                return np.random.choice([ 1, -1])
            else:
                return 0
            
        to_mutate = np.random.binomial(1, self.mutation_frequency, len(vector))
        signs =  map(get_sign, to_mutate)
        mutations = np.array(list(map(lambda x: x*self.mutation_size, signs)))
        return vector + mutations

    def constrain(self, vector):
        def constrainer(x):
            if x < 0:
                return 0
            if x > 1:
                return 1
            return x
        return np.array( list( map(constrainer, vector)))

    def copy(self):
        '''
        Return a new EvolvableVector with mutated parameters.
        '''
        new_vector = self.mutate(self.vector.copy())
        new_vector = self.constrain(new_vector)
        return EvolvableVector(new_vector,
                               self.mutation_size,
                               self.mutation_frequency)
    



class EvolvableVectorConstrained(EvolvableVector):
    """bounding_line is an array of shape (2,X)
    """
    def __init__(self,
                 vector,
                 mutation_size,
                 mutation_frequency,
                 bounding_line):

        EvolvableVector.__init__(self, vector,  mutation_size, mutation_frequency)
        self.bounding_line = bounding_line

    def constrain(self, vector):
        def constrainer(x):
            if x < 0:
                return 0
            if x > 1:
                return 1
            return x

        def bound_to_set(on_line, the_point):
            # Calculate distance to point
            to_line = ((on_line - the_point)**2).sum(axis=1)
            closest_on_line = on_line[to_line.argmin(),:]
            diff = closest_on_line - the_point
            if (diff > 0).all(): # falls above line
                return the_point
            else:
                return closest_on_line
        
        vector = np.array( list( map(constrainer, vector)))
        return bound_to_set(self.bounding_line, vector)

    def copy(self):
        '''
        Return a new EvolvableVector with mutated parameters.
        '''
        new_vector = self.mutate(self.vector.copy())
        new_vector = self.constrain(new_vector)
        return EvolvableVectorConstrained(new_vector,
                                          self.mutation_size,
                                          self.mutation_frequency,
                                          self.bounding_line)
        
