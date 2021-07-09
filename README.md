# Historian

Historian is a python script meant to simulate the double slit experiment using mildly constrained random walks.


### Working hypothesis

The interference patterns from the double slit experiment are conceptually explainable without needing individual particles to have nonlocal characteristics.  The proposed mechanism is a curious sort of Brownian motion.  

The source of this Brownian motion is taken to be the very vacuum fluctuations that are standard quantum lore.  Here this notion is considered literally, and the vacuum is imagined as a chaotic sea of sorts that buffets the particles about.

A particle therefore does not take a straight trajectory, but is chaotically jostled, and the path the particle takes is a jagged, unpredictable affair.  A pollen grain in water moves in a similarly jagged and chaotic manner.

The twist here is that the Brownian motion has some springiness to it, such that the jag lengths are symmetrically distributed around 1 Planck's constant's worth of action.


### Motivation

Feynman's method of summing over all histories is analogous to a typical probability problem.  Figure out all the different paths that can be taken and where they end up.  Count these.  The result is a histogram showing the relative probability of a particle ending up in a particular place.

However, Feynman's method does not simply add the counts together. Instead he weighs each path by the path's "phase".  This detail causes counter intuitive predictions that are nevertheless born out by experiment to astounding accuracy.


#### But what is phase?

Feynman calculates the phase for a particular path in the following way:

```buildoutcfg
phi = exp(-i (2 pi) A / h
```

where A is the path's action, and h is Planck's constant (also in units of action).

Any path where A is precisely an integer number of units h will have a phase of 0 degrees.  Likewise, any path where A is halfway between integer values will have a phase of 180 degrees.  The path's phase, then, is a kind of indicator for the degree to which the path's action is an integer number of Planck units.

Consider a particle with mass m and constant velocity v, traveling a distant d.   It spends an amount of time t with kinetic energy K.  Absent any forces, the path's action is:

```buildoutcfg
A = K t = (m v^2 / 2) (d / v) = m v d / 2 
```

1 Planck's constant worth of action is equivalent to a certain distance:

```buildoutcfg
d = 2 A / mv = 2 h / mv
```

Any path made of precisely 1 unit jags will have a phase of 0.  What's more, any path made of many jags that on average are 1 unit jags will also have a phase of 0.  Paths that deviate from this average will have a phase that deviates from 0. 

Imagine somehow being able to account for every single random walk path the particle could take from the source to every spot on the detector.  Some paths would have many jags that deviate from 1.  Some paths would have few jags that deviate from 1.  If jags of 1 are most likely, then phase 0 paths are also most likely.  The particles would have a tendency to take phase 0 paths, and hit the detector in spots where phase 0 paths congregate.

This is the conceptual reason proposed here behind Feynman's phase based weighting scheme in what is an otherwise standard probability problem.  


#### Simulation results so far

This hypothesis is readily testable by simulating random walks.  The file million_walks.png is the result of one million random walks that successfully hit the detector screen.  It is not an interference pattern.  However, it is not not an interference pattern either.  There is indeed a slight but unmistakable concentration between the slits, and hints of nodes to either side. 


#### Simulation conditions

Each particle begins its history at a point source, and travels one step at a time.  The step is calculated as follows:

- The angle is taken at random.

- The length of the step is sampled from the following distribution:

```buildoutcfg
p = 2 cos^2(pi x)
```
- A plot of this distribution can be found in distribution.png.  It is centered at 1 and tapers to zero at 0.5 and 1.5.

- The integral of this distribution on the interval [0.5, 1.5] is also shown.  It is effectively a quantile function:

```buildoutcfg
q = - 1 / 2 + x + sin(2 pi x) / (2 pi)
```

- If it were solvable, the inverse of this quantile function would be used to map a random number from a uniform distribution to a jag length from the given probability distribution.  

- Instead, the Newton Raphson method is used to approximate the solution, and a table with these values is stored at a resolution of 1000 points.  A random integer is picked between 0 and 1000.  Dividing each by 1000 yields a quantile, and the length x at that quantile is taken as the jag length.

- The file verification.png shows the histogram of 100,000 random numbers processed in this way, reproducing the probability distribution.

- The particle is terminated if a jag takes the particle across a wall or divider boundary, unless it crosses the divider at the location of the slits.

- At the detector screen, the final step of the particle is unlikely to line up exactly with the detector.  Therefore, the final step is chopped off at the detector boundary.  The length of this step is plugged into the probability distribution to get a weighting factor, and this is added to the growing histogram. 


#### Loose ends

The model for this simulation is unsatisfactory in the following ways:

- Obviously, so far the simulation has only produced a subtle interference-ish pattern.  Satisfactory results would need to produce much deeper nodes.  However, there are several parameters at play.  The location of the slits, the distance to the detector, the gap between the slits, these are all free parameters.  Increasing the distance to the back wall has already enhanced the effect compared to smaller dimensions.  The downside is, of course, longer simulation times.

- Feynman's integral assumes a flat detector, and only integrates to this distance.  However, there is no guarantee that a random walk path will land precisely at the detector.  Essentially, only paths that exactly hit the detector should really count, as otherwise they are outside of Feynman's histories.  In order to mimic this effect, the final jag is cut where it meets the detector.  The length of this cut determines the weighting for the history as a whole, as described above.  This is pretty ad hoc, and it will here be noted that without weighting the histories in this way, a very ordinary gaussian distribution appears in the histogram.  

- However, this still serves to highlight that random walks produced in this way lead to a phase based distribution on the detector.  In other words, there are indeed locations where phase 0 paths congregate.

- The cosine squared distribution is convenient, as it is has an analytic integral and has zeros at well defined locations.  This is also a natural choice, as wave equations produce sine waves that then get squared for probability distributions.  However, ideally the pool of possible histories should be the same as Feynman's pool.  Because no jag is less than 1/2 unit, not every Feynman path is represented.  Originally the simulation was meant to use a gaussian distribution around 1 instead, but it was less convenient to work worth.  The most appropriate distribution to use is an open question at this point.

- Feynman's method treats all phases equally.  Really it is the difference in phase between histories that matters.  The model presented here, however, treats phase 0 and phase 180 paths very differently.  This is a noted conceptual mismatch.  For instance, the distance from the source to the detector is an integer number of units, and might be biasing the results.  A successful model should produce a center concentration no matter the distance to the detector.  This is yet to be tested.  

- There is also a lack of direction in this model.  Since each step is random, there is no pervading velocity.  It effectively is a diffusion based model.  There must be some contribution to the pervading direction from particles's velocity. If each random jag is an additive component to a path in a particular direction, instead of the full story, this might resolve the other issues above as well.  I suspect the path of minimum action provides the arbitrary additive constant needed to make the model phase symmetric.  There are plans to test this in the next simulation.


### But what about tunneling?

A successful conceptual model should be able to shed light on the more counter intuitive quantum phenomena.  For example, what about tunneling?

In light of the model presented here, tunneling is simply a matter of Brownian motion kicking the particle beyond where it could go of its own momentum.  Though tunneling is often described as tunneling through a barrier it cannot get over, this is just a metaphor.  There is no barrier to go "over", it is always a matter of going "through".  There is a region of space with higher potential energy that graphically looks like a vertical barrier when plotted with energy on the vertical axis. But physically, the process of going over an energy barrier is really going through some region of space with increased resistance.  Tunneling need not represent some exotic path, only the same path with extra kick behind it.


#### But what about nodes?

Nodes are places with zero predicted particle hits.  In the double slit experiment, opening the second slit seems to place nodes where there were none with one split open.  In other words, locations that were perfectly valid with one slit open are suddenly avoided, even though that original slit is still unblocked.

But a node is infinitesimal in width.  Nodes or no, there are no locations of finite width with zero probability for a particle to hit.  What opening the second slit does, is allow for many new possible paths.  These new paths are not necessarily symmetrically distributed as far as where they lead.  A somewhat likely detector location when one slit is open might not be so likely anymore if there are many new available paths and most of them don't go there.


#### But what about entanglement?

Entanglement is often portrayed in a one to one light, as in every particle that is a member of an entangled pair instantly knows about the other once measurement takes place.  However, to my knowledge, the experimental and theoretical predictions are much more modest, merely a better than chance correlation.  

The conceptual model proposed here allows for enhanced correlations due to the induced clustering around integer units of action.  Two photons created from one event, for instance, are phase correlated.  They have their staring points at the same location, and on average both will move in 1 h units.  The assumption is that other properties would also be affected in 1 h units, and so there would be a tendency for their properties to be correlated.


#### But what about the Uncertainty Principle?

In this model, the uncertainty principle is not so abstract, but a direct consequence of constant jostling.  A particle has an exact location and momentum at any time, but part of this momentum is due to the extra, unknown input from Brownian motion.  Several measurements must be taken to pin down the actual momentum.  But of course, this blows the position measurements.  


#### But what about the measurement paradox?

In this model, there is no measurement paradox, because the particle has a definite though chaotic history.  The wave function has no physical reality, and only represents a histogram.  Wave function collapse is no more physical then rolling a pair of dice while considering all possible outcomes.  Before the dice land, the histogram is the best basis for predicting the results.  But physically it is still particular dice in particular orientations the entire time.  The results appear to be guided by the histogram because after many rolls they will indeed do so.  But nobody would argue that before the dice landed they were just a histogram.  


#### But what about the observer effect?

Of course, the double slit experiment encompasses more than just interference patterns.  The interference goes away upon placing a detector at one slit.

But no detection occurs without interaction, likely with some exchange of momentum.  The particle going through the detected slit is no longer traveling at precisely the same velocity it was before.  Therefore, 1 unit of action is now a different distance, smearing the pattern.  This is why there are weak measurements that disturb the interference pattern less and strong measurements that disturb the pattern more.  Upcoming simulations are planned to test this.


#### References

- QED: The Strange Theory of Matter and Light.  Richard Feynman.
- Quantum Chemistry, Donald A McQuarrie.
- Introduction to Quantum Mechanics.  David J. Griffiths.


#### Thank you!