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
* ability to kill $NIC_1 (PT)$ 
* immune cells constitently do interaction that are either immune-normal or immune-tumor

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

