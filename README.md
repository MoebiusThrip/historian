# Historian

Historian is a python script to simulate the double slit experiment using mildly constrained random walks.

### Motivation

Feynman's method of summing over all possible particle histories is analogous to a typical frequentist probability problem.  Figure out all the different paths that can be taken and where they end up.  Count these.  The result is a histogram showing the relative probability of a particle ending up in a particular place.

However, Feynman's method does not simply add the counts together, instead weighting each path by the path's phase.  This detail causes the counterintuitive, though demonstrably accurate behaviors in the double slit experiment.

#### But what is phase?

Feynmans's method calculates the phase of the path from the path's action.  It is a periodic function of the action.  As the action increases, the phase cycles.  Every Planck's constant increase in action results in a complete cycle.  Therefore, the phase is an indicator for partial Planck units of action.  

A path with integer Planck units of action has phase of 0, no matter how long the actual path is.  A path halfway between integer Planck units of action has a phase of 180.

Therefore, Feynman's method emphasizes locations on the detector where the majority of the paths are phase 0, i.e., reachable by many paths of integer Planck unit action.  

Feynman's method deemphasizes locations reachable on the detector by mixtures of integer and noninteger Planck units of action.  There have phases that cancel out.

#### Working Hypothesis

Suppose there is some medium that randomly affects the motion of the particle, analogous to Brownian Motion.  Unlike Brownian Motion, the probability of these random effects is centered around a Planck unit of action.  Any particular interaction may result in an increase of path length of more than one unit or less than one unit.  But after many such interactions, the path of the particle is likely very close to an integer number of Planck units.  

Therefore, such a medium would produce mainly phase 0 or close to phase 0 paths.  Places on the detector reachable by many phase 0 paths would have a higher probability of being hit.  Places reachable instead by paths of varying phase would be less likely, because only the phase 0 paths have a high probability.  No particle needs to go through both slits, nor interfere with itself, in order to cause an interference pattern on the detector.  

#### The Experiment

Each electron is emitted from a point source, and travels one step at a time.  The step is calculated as follows:

- The angle is taken at random

- The length of the step is sampled from a cos^2 distribution centered at 1 and tapering to 0 at 0.5 and 1.5.  The integral of this distribution can be used to calculate the quantile at all distances.  Reversing this function enables a random number from a uniform distribution to get mapped to a step length, such that they will create the distribution.

- If the step takes the electron across a wall or divider boundary, the electron is terminated. 

- At the detector screen, the final step of the electron is unlikely to line up exactly with the detector.  Therefore, the final step is chopped off at the detector boundary.  The length of this step is plugged into the cos^2 distribution to get a weighting factor for the path. 

#### Results so far

"million_walks.png" shows the latest 
