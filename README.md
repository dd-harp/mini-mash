# mini-mash
the minimum viable simulation

## paper

* RM-delay: one example
* maybe predator-prey with some sort of delay?

## notes

for disaggregated sim with explicit exposed class to track exact output, store the trajectory as a set and then sort it by state change time. Push output *every* time something happens.

It needs to be possible to somehow assign H2M bites (S to E transitions) to S mosquitoes who won’t survive the rest of that TWICE step.

The current disaggregated.cpp has a problem in that mosquitoes which experience a S to E transition experience too much hazard. The part between btime and t0 is experienced twice. This means the hazard for E’s is right but there are too few S’s.

The exact file has another problem. S to E transitions are only assigned to S mosquitoes that are conditioned to survive the part between btime and t0, so the hazard for S mosquitoes is right but there will be too many E’s.

During the BM module, if a S->E mosy transition fires, then we store it. When we use it to update, we can do a little approximation to get things to work out nicely. Let's say we're now in the mosquito module and we get an S->E transition. "That mosquito" (whichever one of the ensemble it occurred to) has already properly accumulated hazard over the TWICE step, so it's technically inappropriate to sample for survival between [btime,t0). However, what we can do is sample P(surv) = exp(-g*(t0-btime)). If they survived, then set S-=1 and E+=1 at t0, because (this is the approximation) that was a S->E transition that landed on a mosquito which survived to the end of the step, t0. If not, don't update anything (it landed on a mosquito that did not survive to the end of the step).

Exact

H2M
Consider the H2M (mosy infection) event. We can use the ensemble trajectory of IH humans. We need the individual life histories of SV mosquitoes. Given IH, we get a time varying hazard (each time IH recovers or is added) that's identical for each mosquito, but they see a different amount of it. Crucially, we only need to sample for those mosquitoes that survive until the end of the TWICE step.

M2H
Because humans can't die, this one is easy. We construct the per-capita hazard trace from IH for each human, and evaluate it for each SH. We need each SH's trajectory, because I->S is an allowed transition so some newly recovered humans may pop up and be eligible for the hazard in the middle of the TWICE step.
