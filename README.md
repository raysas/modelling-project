# Avascular Tumor Growth

## Summary

### Definition of the system

_This is agent-based model in 2D space, Moore neighborhood is used, and Gompertzian growth dynamics is considered_  


There are 2 types of agents:
* NIC: non-immune cells
* IC: immune cells

A group is a collection of agents of the same **type** (NIC/IC).  
Each have different special modes/states; 

for $NIC$:
* 0: normal cell N (or empty cell $\epsilon$, free spaces are like nutrients)
* 1: proliferative tumor cell (PT)
* 2: non-proliferative tumor cell (NT)
* 3: necrotic tumor cell (Ne)

for $IC$:  

* 0: NK cell | natural killer cell, it dies after collision with a tumor cell
* 1: CTL cell | cytotoxic T lymphocyte, it can proliferate and die after collision with a tumor cell, but is cell specific

$\#\ new born\ IC = f (\#\ successful\ defeats)$  (details later, just not that it depends on the number of successes/failures in killing tumor cells)  
Both modes have the same behavior: 
* _random walk_ (moving randomly in grid by substituting an $NIC_0$)
* interact with $NIC_1 (PT)$: anti-tumor, neutra, pro-tumor 
* immune cells constitently do interaction that are either immune-normal (random walk) or immune-tumor (could kill it, be killed or nothing happens)

**Initialization:**  
Starting with a 2d grid of size $L \times L$ (lattice), each cell is represented by an agent, all free cells (NIC 0). Initialize some proliferative cells in the middle and some immune cells at corner, starting with paramters found in table 3 and setting Nmm value (0.2 fig 6, or let user change it - see figure 7-9).

### $NIC$ rules

> [NOTE]
> each cell has an $age$ and $radius$
> at each step, we check the agent mode of agent we're at (PT, NT..) and apply the rules that are type-specific accordingly:

```
# -- this is for all NIC agents that are tumorous
#    normal cells are not considered, no rules

if current_agent.mode = 1
    # PT cell, apply PT rules (in PT section)
    # ....
else if current_agent.mode = 2
    # NT cell, apply NT rules (in NT section)
    # ....
else if current_agent.mode = 3
    # Ne cell, no change (stays in this mode)

```

#### $NIC_1$ | PT

We have 2 types of PT: 1 is mutant, no impact of environment on its proliferation, and 2 is non-mutant and is dependent on the environment. The impact is defined by the number of normal cells surroding it (in Moore neighborhood), the larger the number, the greater is the probability to divide. Biologically it's considering free space as nutrients.

Proliferation probability is a function of time and space $f(age, radius)$, r_p is the division probability, it's equal to:  

* $p_1 = p_0 ( 1 - \frac{r}{Rmax})$ when it's a mutant PT (type 1), $p_0$ is base probablity for a mutant proliferative tumor cell
* $p_2 = \phi_0 \times \#neighboringN \times (1-\frac{r}{Rmax})$ when it's type 2, $\phi_0$ is base probablity for a non-mutant proliferative tumor cell

$Rmax$ is the max radius of the tumor, when it reaches it, it stops growing (dure to lack of nutrients, it's an avascular system) $\implies$ r_p=0

> [!CAUTION]
> I'm not sure what's $\delta_p$ in the following and how it's relate to r_p


```
# behavior of NIC_1, we're at current_agent where:
#   current_agent.agent = NIC
#   current_agent.mode = 1 (i.e., PT)

if current_agent.age > age_threshold OR current_agent.radius > $\delta_p$:
    # -- cell that is too old or radius > deltap 
    #       => quinescence, no more division
    
    current_agent.mode = 2 # -- becomes Non proliferative

else if there exists neighbouring NIC_0 (N) AND r_p > 0:  
    # --this is to divide the cell: the 2 daughter cells will be 
    #   the chosen neighboring N cell and 
    #   the current PT agent itself

    r_1 = probability of choosing that agent
    set chosen_agent (tht is NIC_0)

    chosen_agent.mode = 1 # -- becomes 1st daughter cell of age 0
    chosen_agent.age = 0

    current_agent.age=0 # -- becomes 2nd daughter cell of age 0

else:   #--no neighboring N or no division
    # stays PT (no change)

    current_agent.age ++

fi
```

#### $NIC_2$ | NT

an NT cell can become proliferative if it's back in a particular position (?) (radius < delta p). But on the other hand can die if this radius becomes too large (> dela p + delta n)

```
# behavior of NIC_2, we're at current_agent where:
#   current_agent.agent = NIC
#   current_agent.mode = 2 (i.e., NT)

if current_agent.radius < delta_p:
    current_agent.mode = 1 # -- becomes proliferative

```


### $IC$ rules

Immune cells move towards the center of the tumor (Ne debris signals). It starts from one corner of the lattice and moves towards the center.

This movement, _random walk_, is associated with a probability r_walk.
$$r\_walk = k \times [ \frac{nT}{ncell^2}, \frac{nPT}{nT}]$$  

> [!CAUTION]
> not sure how is this calculated

At the start, $\frac{nPT}{nT} \gt \gt \frac{nT}{ncell^2} \iff$ tumor cells are more proliferative, movement starts biased. Later,  $\frac{nT}{ncell^2}$ plays a role in the mvt, it's unbiased, random.

There are 3 possible behaviors when IC meets an NIC of mode=1 (only interacts with PT):

1. it kills the PT cell by replacing it (NIC 1 disappears), here it depends whether the IC is of mode 1 or 2
    * if it's NL, it does the same for all neighboring PT cells
    * else it's CTL (2), it chooses one randomly and does the transition
2. no change, keeps earching
3. it dies, becoming N or $\epsilon$ (agent type NIC, mode 0), and the cancer cell survives.

if there is no neighboring PT, it proceeds with the random walk (searchs for free space to move to it in the next iteration, if there arent any it stays)

recruitement, as we said is a function of success/failure difference in killing tumor cells:
$$\#newbornIC=(v-f) \frac{nPT}{nT}$$

if >0, there is recruitement, which means, a random free space searched becomes IC, and nb cells searching = # newborn IC

> [!CAUTION]
cells only signals other cells that have physical contact with?

### Notes

Some useful functions to implement (get the dynamic param in table 2 or even compute probas)

On an agent level:
- is_immune(agent) -> returns true if current.agent is IC
useful to compute $nI$ and $nI_1$
- is_tumor(agent) -> returns true if agent is NIC and mode is either 1, 2 or 3 (According to figure 4, Ne + NT + PT =total of tumor cells (T) so $nT=Ne+NT+PT$)

at each time step, on the grid:
- avg_radius
- growth_fraction (GF) = $\frac{nPT}{nT}$ (nT here signifies the whole tumor)
- necrotic_fraction (NF) = $\frac{nNe}{nT}$
- tumor_radius (Rt)= avg radiu sof all tumor cells
- necrotic_radius (Rne)= avg radius of all necrotic cells 
- Nmm_t(all) fraction of non mutant to mutant cells

Parameters:
- Nmm should be set at teh beginning of the simulation, that is, tha ration of non-mutant to mutant cells. maybe make the user choose it from a list of values (bar)

The first 9 figures are values of parameters on a model of only NIC agents (still no IC introduced), maybe better to start like this and make sure the same behavior of groth holds, before starting with interactions.

The figrues show:
- 4: number of cells over time (3 types of tumor cells) | not part of the results
- 5: average radius change over time (overall-tumor and necrotic radii) | growth and necrotic fractions over time
- 6: number of mutant and non-mutant cells over time (starting Nmm=0.2)
- 7: trying different values of Nmm (`range(0, 1, step=0.1)`) and plotting the number of cells over time (3 types of cells, each in a seperate plot)
- 8: on different calues of Nmm, GF and NF over time
- 9: on different values of Nmm (initial values), mutant/non-mutant evolution over time (inverse Nmm_t maybe). I think the y label in the figure is wrong

