# disaggregated simulation algorithm

Recall that the "big picture" algorithm look like this:

```
tmax <- t0 + dt
while (tnow < tmax) {

    simulate_mosquitio_pop(mosquito_pop, tnow, dt)
    simulate_human_pop(human_pop, tnow, dt)

    bloodmeal(human_pop, mosquito_pop, tnow, dt)
    
    tnow += dt
}
```

## mosquito simulation on first time step

Consider the problem of starting with simulation of the mosquito population `simulate_mosquitio_pop`, at the initial time point. Then there is no input of bites yet, so the algorithm can directly begin to simulation $S_{V},E_{V},I_{V}$. Let biting input be stored in a queue data structure called `H2M_bloodmeal` (human to mosquito transmission).

Because there is no input of bites which would cause an $S\rightarrow E$ event, we can directly simulate the events in the system. 

If $E_{V}$ is nonzero, we queue $E\rightarrow I$ events for each exposed mosquito, call the queue data structure which stores these events `E2I_mosquito`. This is done "before simulation starts". The time that each of these mosquitos will become infected is a realization of the random variable $Unif(-EIP,0)+EIP$. A draw is taken for each, and it is added to `E2I_mosquito`. Values are drawn from this random variable because if a mosquito was incubating at time $0$, it could not have become infected at a time less than $-EIP$, otherwise it would have been infectious at time $0$. Therefore, under the assumption that the model is (except for the deterministic delays), a set of coupled Poisson processes (because the infection time would have arisen from one), the actual infection time is uniformly distributed along the interval $(-EIP,0]$. Therefore the $E\rightarrow I$ transition will be that initial infection time plus $EIP$.

After this sampling of initial conditions, the mosquito population can be simulated along the time interval from $[t_{0},t_{0} + \Delta t)$. Henceforth let $t_{0} + \Delta t = t_{1}$, so that generally the $i^{th}$ iteration is denoted by its right open endpoint $t_{i}$. Deaths from each compartment occur at per-capita intensity $g$, so for example the Poisson process accounting for deaths from $S_{V}$ has intensity $g S_{v}$. Emergence is even simpler, as we assume it occurs with state-independent intensity $\lambda$. If the next event is an $E\rightarrow I$ transition, we must update slightly more carefully, by removing it from the head of the queue `E2I_mosquito` during the state update. If the event is a death from $E_{V}$, then a random element from `E2I_mosquito` is removed (note that the size of the queue is always equal to the value of $E_{V}$). The whole system can be simulated by Anderson's "modified next reaction method" (MNRM), a technique to simulate Markovian or non-Markovian continuous time discrete event stochastic systems. The mosquito state is simulated using the MNRM until the next event time would exceed $t_{1}$, the end of the TWICE interval. 

During this whole time we must also accumulate the piecewise constant (cadlag) trajectory describing the changes in $\tilde{S}_{V}$ over this initial time interval $[t_{0},t_{1})$. Note also that mosquitos in this trajectory represent mosquitos which may be either in $S_{V}$ or $E_{V}$ (hence the tilde notation), because the bloodmeal events have not been sampled yet, though we are always guaranteed that at $t_{0}$ all the mosquitos are in $\tilde{S}_{V}$. Also note that emergence events which increase $\tilde{S}_{V}$ by 1 will increase the population of susceptible mosquitos, although deaths, which decrease $\tilde{S}_{V}$ by 1 may either occur to a susceptible mosquito or an infected mosquito. This will be relevant when sampling bloodmeal events later. 

The mosquitos must also output the trajectory of $I_{V}$ over this interval, needed for the computation of bloodmeals causing infection in humans.

## human simulation on first time step

Because we do not model deaths or births in the human component, a description of the algorithm for `simulate_human_pop` is even simpler. Again, for non-zero $E_{H}$, we do a similar process as to the initial mosquito incubating individuals to sample the initial state. The difference is that the random variable which is sampled and added to `E2I_human` is now described by $Unif(-LEP,0)+LEP$, following the same argument.

The only Poisson process in the human model is that governing the $I\rightarrow S$ transition, firing at intensity $r I_{H}$. It competes with the deterministic events queued in `E2I_human`. The same update procedure based on the MNRM used for the mosquitos can be used to simulate the humans.

Like the mosquitos, the humans must output the trajectory of $\tilde{S}_{H}$ over this interval. They must also output the trajectory of $X$ (prevalence), needed for the computation of bloodmeals causing infection in mosquitos.

## bloodmeal algorithm on the first time step

### infection in mosquitos $S_{V}\rightarrow E_{V}$ on first time step

We now need to introduce some notation. In the following notation, the subscript $i$ will denote the right hand open endpoint for a piecewise constant interval in the input trajectory of susceptible mosquitos. Let $q_{0}=t_{0}$, and $q_{1}$ is the right endpoint of the first interval. Let $q_{n}$ denote the right hand open endpoint of the final interval, which will equal $t_{1}$. Now, let $s_{i}$ denote the value of the state trajectory $\tilde{S}_{V}$ in this interval. Let $\tilde{r}_{i}$ denote the _putative risk set_ over this interval, that is the number of susceptible mosquitos at the left hand closed endpoint $q_{i-1}$. Let $r_{i}$ denote the _risk set_ at the end of this interval, that is, the number of susceptible mosquitos at right hand open endpoint $q_{i}$.

On the first step $i=1$ the putative risk set is $\tilde{r}_{1} = s_{1} = \tilde{S}_{V}(t), t\in [q_{0},q_{1})$. From the human module we know the value of prevalence $X$ over that interval as well. Each mosquito in the putative risk set has probability $p_{i}$ of taking a bloodmeal that results in infection during that interval, calaculated as:

$$
p_{i} = 1 - \exp\left( -\int_{q_{0}}^{q_{1}} acX(s) ds \right)
$$

Then the number of infections in mosquitos that become infected in the first interval is a draw from the random variable:

$$
n_{i} \sim Binomial(\tilde{r}_{i},p_{i})
$$

Then the risk set (number of mosquitos still susceptible) at the end of the interval ($q_{1}$) is $r_{i} = \tilde{r}_{i} - n_{i}$. Due to the assumption of the model as a set of competing Poisson processes, each of the $n_{i}$ infection events occurs at a time uniformly distributed on the interval, that is $b_{i} \sim Unif(q_{0},q_{1}$$. Each of these bites must be added to the `H2M_bloodmeal` queue.

### infection in mosquitos $S_{V}\rightarrow E_{V}$ on subsequent time steps

When $i > 1$, the algorithm is similar, but must account for the changing size of the risk set due to demographic events in the mosquito population, which are the jumps in $\tilde{S}_{V}$. At the moment in time $q_{i-1}$ the trajectory jumps from $s_{i-1}$ to $s_{i}$, the difference is given by $s_{i} - s_{i-1}$. Because the model is a set of simple (not compound) counting processes, the difference can only take values of $+1$ or $-1$. 

Let's see how this influences the size of the risk set. The size of the risk set _right before_ the jump is $r_{i-1}$. If the jump took a value of $+1$ at $q_{i-1}$, then we know the risk set grew by $1$, because we assume all mosquitos emerge susceptible, and $\tilde{r}_{i} = r_{i-1} + 1$. If the jump took a value of $-1$, then we have to sample to determine if the death occurred among the set of infected or susceptible mosquitoes at $q_{i-1}$. The probability the death occurred to a susceptible mosquito (i.e. a mosquito in the risk set) is $r_{i-1}/s_{i-1}$, and the probability it occurred to an infected mosquito is $1-r_{i-1}/s_{i-1}$. Sample a Bernoulli random variate with those probabilities, if the death occurred to a susceptible mosquito, the risk set decreases by $1$, that is $\tilde{r}_{i} = r_{i-1} - 1$. If the death occurred to an infected mosquito, the risk set does not decrease, that is $\tilde{r}_{i} = r_{i-1}$. Note that we must now "update the history", however we are tracking the overall trajectory of the model, and assign that death concretely to the susceptible or infected population as appropriate (i.e. $S_{V}$ or $E_{V}$).

Now that we know the putative risk set, we can again calculate the per-capita probability of infection $p_{i}$, and sample $n_{i}$.


























We show the algorithm below in pseudocode. Let `s` equal $\tilde{S}_{V}$, which can be linearly indexed by interval.

```
s_prev <- s[1]
r_prev <- 0
i <- 1
while (i <= n)

  s <- s[i]
  
  if s - s_prev == -1
    if rand() < r_prev/s_prev # death occured among susceptible
      rbar <- r_prev - 1
    else
      rbar <- r_prev
      update_history() # reassign this death event to decrease the E_{V} population
    end
  else if s - s_prev == 1
    rbar <- r_prev + 1
  else if s - s_prev == 0
    if i != 1
      error("all interior intervals must correspond to an emergence or death event")
    else
      rbar <- s
    end
  else 
    error("jumps larger than one are not allowed")
  end
  
  r_prev <- r
end
```

### infection in humans $S_{H}\rightarrow E_{H}$
