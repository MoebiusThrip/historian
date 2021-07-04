# Historian

Historian is a python script to simulate the double slit experiment using mildly constrained random walks.


### Working hypothesis

The interference patterns from the double slit experiment are conceptually explainable without needing individual particles to have nonlocal characteristics.  The proposed mechanism is a curious sort of Brownian motion.  

The source of this Brownian motion is taken to be the very vacuum fluctuations that are standard quantum lore.  Here this notion is considered literally, and the vacuum is imagined as a chaotic sea of sorts that buffets the particles about.

A particle therefore does not take a straight trajectory, but is chaotically jostled, and the path the particle takes is a jagged affair.  A pollen grain in water moves in a similarly jagged and chaotic manner.  

The twist here is an unusual pattern in the Brownian motion, such that the distribution of jags in the particle's path is centered symmetrically around 1 Planck unit of action, h.


### Motivation

Feynman's method of summing over all possible particle histories is analogous to a typical probability problem.  Figure out all the different paths that can be taken and where they end up.  Count these.  The result is a histogram showing the relative probability of a particle ending up in a particular place.

However, Feynman's method does not simply add the counts together. Instead he weighs each path by the path's phase.  This detail causes counter intuitive predictions that are nevertheless born out by experiment to astounding accuracy.


#### But what is phase?

Feynman calculates the phase for a particular path in the following way:

```buildoutcfg
phi = exp(-i (2 pi) A / h
```

where A is the path's action, and h is Planck's constant (itself also in units of action).

Any path where A is precisely an integer number of units h will have a phase of 0 degrees.  Likewise, any path where A is halfway between integer values will have a phase of 180 degrees.  The path's phase, then, is a kind of indicator for the degree to which the path's action is an integer number of Planck units.

Any path made of precisely 1 unit jags will have a phase of 0.  What's more, any path made of jags that on average are 1 unit will also have a phase of 0.  Paths with half integer action values would be much less likely.

Imagine somehow being able to account for every single path the particle could take.  The phase 0 paths would be most likely, and especially more likely than the phase 180 paths.  Imagine looking at a spot on the detector, and being able to see all the possible paths that lead from the source to that spot.  If a high proportion of those paths are phase 0 paths, it is a likely place for the particle to go.  If a low proportion of those paths are phase 0 paths, it is a less likely place for the particle to go.

This is the proposed conceptual reason for Feynman's phase based weighting scheme in what is an otherwise standard probability problem.  


#### Simulation results so far

This hypothesis is readily testable by simulating these random walks and counting where they land.  The file "million_walks.png" is the result of one million successful random walks.  It is not an interference pattern.  However, it is not not an interference pattern either.  There is indeed a slight but unmistakable concentration between the slits, and the barest hints of nodes to either side. 


#### The experiment

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

- Instead, the Newton Rhapson method is used to approximate the solution, and a table with these values is stored at a resolution of 1000 points.  Thus, a random integer is picked between 0 and 1000.  Dividing each by 1000 yields a quantile, and the length x at that quantile is taken as the jag length.

- The file "verification.png" shows the histogram of 100,000 random numbers processed in this way, reproducing the cosine squared distribution.

- If a jag takes a particle across a wall or divider boundary, the particle is terminated.  

- At the detector screen, the final step of the particle is unlikely to line up exactly with the detector.  Therefore, the final step is chopped off at the detector boundary.  The length of this step is plugged into the cos^2 distribution to get a weighting factor for the path. 


#### Loose ends

The model for this experiment is unsatisfactory in the following ways:

- Obviously, the simulation only produces a subtle interference-ish pattern.  Satisfactory results would need to produce much deeper nodes.  However, there are several parameters at play.  The location of the slits, the distance to the detector, the gap between the slits, these are all free parameters.  Increasing the distance to the back wall has already enhanced the effect compared to smaller setups.  The downside is, of course, longer simulation times.

- Feynman's integral assumes a flat detector.  However, there is no gaurantee that a random walk path will land precisely at the detector.  Essentially, only paths that exactly hit the detector should really count.  In order to mimic this effect, the final jag is cut when it meets the detector.  The length of this cut determines the weighting for the path as a whole.  Thus, paths that most closely land on the detector are given a higher weight than those 1/2 h away.  For convenience, the same cosine squared distribution is used for this as well, but with no justification.

- This is somewhat ad hoc, and it will here be noted that without correcting the paths in this way (any particle that crosses the detector counts equally), a very ordinary gaussian distribution appears.  However, this still serves to highlight the fact that random walks produced in this way lead to a phase based distribution on the detector.  In other words, there are indeed locations that favor phase 0 paths, despite every jag being a random event.

- The cosine squared distribution is convenient, as it is analytically integratable and has zeros at well defined locations.  This is also a natural choice, as wave equations produce sine waves that then get squared for probability distributions.  However, ideally the pool of possible paths should mimic the pool of Feynman's paths.  Because no jag is less than 1/2 h, it is uncertain that every Feynman path is represented.  i.e., there should be a chance for all paths, no matter how unlikely.  Other distributions could be tried.

- Feynman's method treats all phases equally.  Really it is the difference in phase between paths that matters.  The model presented here, however, treats phase 0 and phase 180 paths very differently.  This is a noted conceptual mismatch.

- Also, the model lacks any notion of direction.  All particles are assumed to have the same momentum, allowing a mapping of jag length to action.  However, there is no way for a particle to move in a direction over time, as at each jag the angle is random.  

- A more sophisticated approach could grant the particles overall momentum in some direction, with the jags only responsible for deviations from this path.  In this way, the sum total of deviations would be 0, though the path may have an overall nonzereo phase.  Here it will be noted that in bound states, it is only in the classical portion of the wavefunction that the phase changes.  Far from the nucleus, the kinetic energy of the particle becomes negative, and there is no phase change.  Despite phase being a quantum mechanical concept, changes in phase are strictly classical.  


### But what about tunneling?

A successful conceptual model should be able to shed light on the more bizarre of quantum phenomena.  In light of the model presented here, tunneling is simply a matter of Brownian motion kicking the particle beyond where it could go of its own momentum.  Though tunneling is often described as tunneling through a barrier it cannot get over, this is just a metaphor.  There is no barrier to go over.  There is a region of space with higher potential energy that graphically looks like a vertical barrier when plotted with energy on the vertical axis. But physically, the process of going over an energy barrier is really going through some region of space with increased resistance, tunneling or no.


#### But what about nodes?

Nodes are places with zero predicted particle hits.  In the double slit experiment, opening the second slit seems to place nodes where there were none with one split open.  In other words, locations that were perfectly valid with one slit are suddenly avoided, even though that original slit is still unblocked.

But a node is infinitesimal in width.  Nodes or no, there are no locations of finite width with zero probability for a particle to hit.  What opening the second slit does, however, is allow for many other possible paths.  But these new paths are not necessarilily symmetrically distributed as far as where they lead.  A somewhat likely detector location when one slit is open might not be so likely anymore if most of the new possible paths don't go there.


#### But what about entanglement?

Entanglement is often portrayed in a one to one light, as in every particle that is a member of an entangled pair instantly knows about the other once measurement takes place.  However, to my knowledge, the experimental and theoretical predictations are not one to one, but merely better than chance.  The conceptual model proposed here allows for enhanced correlations due to the induced clustering around integer units of action.  Two photons created from one event, for instance, are phase correlated.  They have their staring points at the same location, and on average both will move in 1 h units.  Other properties would also be affected on average in 1 h units, and so there would be a tendency for their properties to be correlated.  But with no communication.  Only because the Brownian motion in the intervening space tends toward 1 h, phase neutral effects.


#### But what about the Uncertainty Principle?

In this model, the uncertainty principle is not so abstract, but a direct consequence of constant jostling.  A particle has an exact location and momentum at any time, but part of this momentum is due to the extra, unknown input from Brownian motion.  Several measurements must be taken to pin down the actual momentum.  But of course, this blows the position measurements.  


#### But what about the measurement paradox?

In this model, there is no measurement paradox, because the particle has a definite though chaotic path.  The wave function has no physical reality, and only represents a histogram.  Wave function collapse is no more physical then rolling a pair of dice while considering all possible outcomes.  Once the dice lands, you no longer need to think about the histogram of possibilities.  


#### But what about the observer effect?

This simulation allows for testing the observer effect as well.  No detection can occur without some kind of interaction.  If this interaction adds a random phase to the particle's path, it would be expected to smear out the interference pattern.  A stronger interaction would be expected to more effectively scatter the phase.  This is why there are both strong and weak measurements that erase the inference pattern to different degrees.


#### References

- Quantum Chemistry, Donald A McQuarrie.
- Introduction to Quantum Mechanics.  David J. Giffiths.
- QED: The Strange Theory of Matter and Light.  Richard Feynman.


#### Thank you!