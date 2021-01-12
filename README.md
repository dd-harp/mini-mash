# Mini-MASH
The minimum Ross-Macdonald (RM) style model for which a disaggregated representation makes sense. The model is **RM-delay**, because E&rarr;I transitions in both the human and mosquito populations only occur after a fixed delay. In this case a disaggregated representation means that different portions of code are responsible for the various components of the model. The components are disjoint subsets of the full state space such that their union gives the full state space of the model (an aggregated representation).

The key point is that feedback loops in the standard compartmental flow diagram representation of the system can be broken across time steps if there is a minimum delay in the feedback loop, such that a time step can be chosen that is less than or equal to that minimum delay. We call such a time step TWICE, the Temporal Window of Indifference to Contingent Events. In this case, an additional component that samples events rather than stores states should be introduced that is responsible for sampling the initiation times of events with delay. For **RM-delay**, it is the bloodmeal.

Causal dependency between mosquito and human populations can only be transmitted via the bloodmeal (when pathogens are exchanged). During a TWICE interval, bloodmeals which cause S&rarr;E events in humans and mosquitoes are sampled. Because by definition a TWICE interval must be less than the LEP (Liver Emergence Period) and EIP (Extrinsic Incubation Period), those "initiating" S&rarr;E events will not result in E&rarr;I events until the next TWICE step. Because only I individuals can affect the other species, given the system state at the start of the TWICE interval, the two populations are conditionally independent until the next one; that is, they can be simulated independently until the start of the next TWICE interval, where they must find out what S&rarr;E bloodmeals happened on the last step.

## Directory structure

* data: output from simulation runs
* deprecated: unmaintained files
* figs: produced figures
* scratch: test files
* scripts: run scripts (prefaced with run) and figure making scripts (prefaced with fig)
* src: simulation source code
