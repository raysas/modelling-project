'''
Replication:
Shirakawa, T., Sato, H., & Ishiguro, S. (2015). 
Construction of living cellular automata using the Physarum plasmodium. 
International Journal of General Systems, 44(3), 292-304. 
https://doi.org/10.1080/03081079.2014.997531
'''

import os,sys
sys.path.append(os.path.abspath('Cellular-Automata-Python'))

from cellularautomata import CountType, GuiCA
from random import random

# -- taken from table 1 p7
state_transition_probabilities = {
    (0, 0): {1: 0.987, 2: 0.976, 3: 0.954, 4: 0.979,0: 1},
    (0, 1): {1: 0.013, 2: 0.024, 3: 0.046, 4: 0.021, 0:0},
    (1, 0): {0: 0.124, 1: 0.041, 2: 0.016, 3: 0.008, 4: 0.002},
    (1, 1): {0: 0.876, 1: 0.959, 2: 0.984, 3: 0.992, 4: 0.998},
}



def physarum(cell, neighbors: list):
    '''
    shirakawa (2015) model 1 implementtaion
    make sure to set von neumann neighborhood
    '''
    state, time= cell
    num_state_one_neighbors = CountType(neighbors, 1)
    # print(f'cell: {cell}, neighbors: {neighbors}, num_state_one_neighbors: {num_state_one_neighbors}')
    match state:
        case 'zero': 
            # -- need to cehck if has at least 1 state one neighbor (table 1 p7)
            proba_random = random()
            if proba_random < state_transition_probabilities[(0, 1)][num_state_one_neighbors]:
                return ('one', None)
            else:
                return ('zero', None)
            
        case 'one':
            proba_random = random()
            if proba_random < state_transition_probabilities[(1, 0)][num_state_one_neighbors]:
                # print(f'probability of 1 to 0: {state_transition_probabilities[(1, 0)][num_state_one_neighbors]} -> going to 0')
                return ('zero', None)
            else:
                return ('one', None)
            
        case _:
            return cell
            

def main():
    cellcolors = {   
        ('zero', None): 'white',
        ('one', None): 'yellow',
        }

    GuiCA(physarum, cellcolors)

if __name__=='__main__':
    main()