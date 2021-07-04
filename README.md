# Historian

Historian is a python script to simulate the double slit experiment using mildly constrained random walks.


### Working Hypothesis

The interference patterns from the double slit experiment are conceptually explainable without needing individual particles to have nonlocal characteristics.  The proposed mechanism is a curious sort of Brownian motion.  

The source of this Brownian motion is taken to be the very vacuum fluctuations that are standard quantum lore.  Here this notion is considered literally, and the vacuum is imagined as a chaotic sea of sorts that buffets the particles about.

A particle therefore does not take a straight trajectory, but is chaotically jostled, and the path the particle takes is a jagged affair.  A pollen grain in water moves in a similarly jagged and chaotic manner.  

The twist here is an unusual pattern in the Brownian motion, such that the distribution of jags in the particle's path is centered symmetrically around 1 Planck unit of action, h.


### Motivation

Feynman's method of summing over all possible particle histories is analogous to a typical probability problem.  Figure out all the different paths that can be taken and where they end up.  Count these.  The result is a histogram showing the relative probability of a particle ending up in a particular place.

However, Feynman's method does not simply add the counts together. Instead he weighs each path by the path's phase.  This detail causes some counter intuitive results that are nevertheless born out by experiment to astounding accuracy.


#### But what is phase?

Feynman calculates the phase for a particular path in the following way:

```buildoutcfg
phi = exp(-i (2 pi) A / h
```

where A is the path's action, and h is Planck's constant (itself also in units of action).

Any path where A is precisely an integer number of units h will have a phase of 0 degrees.  Likewise, any path where A is halfway between integer values will have a phase of 180 degrees.  The path's phase, then, is a kind of counter for the degree to which the path's action is an integer number of Planck units.

Any path made of precisely 1 unit jags will have a phase of 0.  What's more, any path made of jags that on average are 1 unit will also have a phase of 0.  Paths with half integer action values would be much less likely.

Imagine somehow being able to account for every single path the particle could take.  The phase 0 paths would be most likely, and especially more likely than the phase 180 paths.  Imagine looking at a spot on the detector, and being able to see all the possible paths that lead to that spot.  If a high proportion of those paths are phase 0 paths, it is a likely place for the particle to go.  If a low proportion of those paths are phase 0 paths, it is a less likely place for the particle to go.

This is the proposed conceptual reason for Feynman's phase based weighting scheme in what is an otherwise standard probability problem.  


#### Simulation Results So Far

This hypothesis is readily testable by simulating these random walks and counting where they land.  The file "million_walks.png" is the result of one million successful random walks.  It is not an interference pattern.  However, it is not not an interference pattern either, as there is indeed a slight but unmistakable concentration between the slits, and the barest hints of nodes to either side. 


#### The Experiment

Each particle is emitted from a point source, and travels one step at a time.  The step is calculated as follows:

- The angle is taken at random

- The length of the step is sampled from the following distribution:

```buildoutcfg
p = 2 * cos^2(pi * x)
```
- A plot of this distribution can be found in distribution.png.  It is centered at 1 and tapers to zero at 0.5 and 1.5.

- The integral of this distribution on the interval [0.5, 1.5] is also shown.  The integral is:

```buildoutcfg
q = x - (1/2) + sin(2 * pi * x) / (2 * pi)
```

- If it were solvable, the inverse of this integral would be used to map a random number from a uniform distribution to a jag length from the squared cosine distribution given above.  

- Instead, the Newton Rhapson method is used to approximate the solution, and a table with these values is stored at a resolution of 1000 points.  Thus, a random integer is picked between 0 and 1000.  Dividing by 1000 yields a quantile, and the length x at that quantile is taken as the jag length.

- The file "verification.png" shows the histogram of 100,000 random numbers processed in this way, reproducing the cosine squared distribution.

- If a jag takes a particle across a wall or divider boundary, the particle is terminated.  

- At the detector screen, the final step of the particle is unlikely to line up exactly with the detector.  Therefore, the final step is chopped off at the detector boundary.  The length of this step is plugged into the cos^2 distribution to get a weighting factor for the path. 


#### Unsolved Issues



#### But What About Nodes?

#### But What About Entanglement?

#### But What About the Uncertainty Principle?

#### But What About the Measurement Paradox?

#### But What About the Observer Effect?


#### References

