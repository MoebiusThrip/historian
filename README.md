# Historian

Historian is a python script meant to simulate the double slit experiment using guided random walks.


### Working hypothesis

The interference patterns from the double slit experiment are conceptually explainable without needing individual particles to have nonlocal characteristics.  The proposed mechanism is a curious sort of Brownian motion.  

The source of this Brownian motion is taken to be the very vacuum fluctuations that are standard quantum lore.  Here this notion is considered literally, and the vacuum is imagined as a chaotic sea of sorts that buffets the particles about.

A particle therefore does not take a straight trajectory, but is chaotically jostled, and the path the particle takes is a zigzagged, unpredictable affair.  A pollen grain in water moves in a similarly zigzagged and chaotic manner.

The twist here is that the Brownian motion has some springiness to it, such that very small zigzag lengths are not as likely as zigzags some distance away.


### Motivation

Feynman's method of summing over all histories is analogous to a typical probability problem.  Figure out all the different paths that can be taken and where they end up.  Count these.  The result is a histogram showing the relative probability of a particle ending up in a particular place.

However, Feynman's method does not simply add the counts together. Instead it weighs each path by the path's "phase".  This detail causes counter intuitive predictions that are nevertheless born out by experiment to astounding accuracy.


### But what is phase?

Feynman calculates the phase for a particular path in the following way:

```buildoutcfg
phi = exp(-i (2 pi) A / h
```

where A is the path's action, and h is Planck's constant (also in units of action).

Despite its mysterious connotations the imaginary i in the exponent simply describes a periodic function.  Any path where A is precisely an integer number of units h will have a phase of 0 degrees.  Likewise, any path where A is halfway between integer values will have a phase of 180 degrees.  The path's phase, then, represents the remainder after accounting for all full Planck's constants' worth of action.

Consider a particle with mass m and constant velocity v, traveling a distance d.   It spends an amount of time t with kinetic energy K.  Absent any forces, the path's action is:

```buildoutcfg
A = K t = (m v^2 / 2) (d / v) = m v d / 2 
```

1 plancksworth of action at a particular velocity is equivalent to a certain distance:

```buildoutcfg
1 pw = 2 A / m v = 2 h / m v
```

Consider the path of a particle flying towards a screen that measures 100.1 pw.  Consider a slightly different path that measures 101.1 pw.  Both these paths would have a phase of 36 degrees.  They are "in phase".

Consider a third path that measures 61543009.1 pw.  Also 36 degrees and also in phase.  In fact, adding any integer number of plankworths retains the in phase relationship.  Adding any number of zigzags of precisely one planksworth keeps the paths in phase.

What's more, adding any number of zigzags, whose average length is one planksworth will also tend to maintain the in phase relationship, especially once the number of zigzags is large.

Imagine you are approaching a row of empty chairs, seeking one in which to sit.  Now imagine your feet are strapped into a snowboard.  Now imagine you are drunk, and really need to sit down.  It may not be obvious from far away, but no matter how you stumble about there will be certain chairs that are easier to get to.  A chair next to an easy chair will likely not be an easy chair, rather the easy chairs will be spaced at intervals in accordance with the size of your snowboard. 

This is the proposed connection with the Feynman weighting scheme.  There are some places on the detection screen that are simply easier to get to than others.  Paths that are in phase are compatible with each other, because unit zigzags are easy.  Feynman's method picks out the locations on the detector most easily reachable by random zigzags of unit steps.


### Simulation results so far

This hypothesis is readily testable by simulating these random walks.  The file million_walks.png is the latest result of one million random walks that successfully hit the detector screen.  Certainly it is not an obvious interference pattern.  And yet I daresay there is a hint.


### Simulation conditions

Each particle begins its history at a point source, and travels one zigzag at a time.  The zigzag is calculated in two parts:

The first step is a zig.  It represents classical motion towards the detection screen along a straight path.  The length of this step is taken to be random.  In practical terms, this means picking a random number from a uniform distribution up to some maximum length, say 0.5 planckworths.

The second step is a zag.  It represents the quantum component, and is calculated in a more complicated way.

The angle is taken at random.

The length of the step is sampled from the following distribution:

```buildoutcfg
p = 2 cos^2(pi x)
```
A plot of this distribution can be found in distribution.png.  It is centered at 1.0 and tapers to zero at 0.5 and 1.5.

The integral of this distribution on the interval [0.5, 1.5] is also shown.  It is effectively a quantile function, in that it calculates the fraction of samples from the above distribution with a length smaller than that given.  Applying the function to 1.0, for instance, results in a quantile of 0.5 because it is at the median of the distribution.

```buildoutcfg
q = - 1 / 2 + x + sin(2 pi x) / (2 pi)
```

If it were solvable, the inverse of this quantile function would be used to map a random number from a uniform distribution to a zag length from the given probability distribution.  

Instead, the Newton Raphson method is used to approximate the solution, and a table with these values is stored at a resolution of 1000 points.  A random integer is picked between 0 and 1000.  Dividing each by 1000 yields a quantile, and the length x at that quantile is taken as the zag length.

The file verification.png shows the histogram of 100,000 random numbers processed in this way, reproducing the probability distribution.

- The particle is terminated if a zig or zag takes the particle across a wall or divider boundary, unless it crosses the divider at the location of the slits.

- At the detector screen, the final step of the particle is unlikely to line up exactly with the detector.  Therefore, the final step is chopped off at the detector boundary, and the detection is weighed based on the likelihood of a zag that length.  A zig is granted full weight no matter how long, as the distribution of zig lengths is taken to be uniform. 


### Loose ends

The model for this simulation is unsatisfactory in the following ways:

- Obviously, so far the simulation has only produced a subtle interference-ish pattern.  Satisfactory results would need to produce much deeper nodes.  However, there are several parameters at play.  The location of the slits, the distance to the detector, the gap between the slits, the frequency of zags, the zag distribution function, the detection weighting function, these are all free parameters.  Currently running simulations are exploring this parameter space.

- Feynman's integral assumes a flat detector, and only integrates to this distance.  However, there is no guarantee that a random walk path will land precisely at the detector screen.  Essentially, only paths that exactly hit the detector should really count, as otherwise they are outside of Feynman's histories.  In order to mimic this effect, the final zag is cut where it meets the detector.  The length of this cut determines the weighting for the history as a whole, as described above.  This is pretty ad hoc, and it will here be noted that without weighting the histories in this way, a very ordinary gaussian distribution appears in the histogram.  

- However, this still serves to highlight that random walks produced in this way lead to a phase based distribution at the detector.  In other words, there are indeed locations where in phase paths congregate.

- The cosine squared distribution is convenient, as it is has an analytic integral and has zeros at well defined locations.  This is also a natural choice, as wave equations produce sine waves that then get squared for probability distributions.  However, ideally the pool of possible histories should be the same as Feynman's pool.  Because no zag is less than 1/2 unit, not every Feynman path is represented.  Originally the simulation was meant to use a gaussian distribution around 1 instead, but it was less convenient to work worth.  The most appropriate distribution to use is an open question at this point.


### Further Speculations

A successful conceptual model should also be able to shed light on the more counter intuitive quantum phenomena.  For example, what about tunneling?


#### What about tunneling?

In light of the model presented here, tunneling is simply a matter of Brownian motion kicking the particle beyond where it could go of its own momentum.  Though tunneling is often depicted as tunneling through a barrier it cannot get over, this is just a metaphor.  Tunneling does not represent an alternative path, but rather further along the same path.  There is no barrier to go "over", it is always a matter of going "through".  There is a region of space with higher potential energy that graphically looks like a vertical barrier when plotted with energy on the vertical axis. But physically, the process of going over an energy barrier is really going through some region of space with increased resistance.  Tunneling need not represent some exotic path, only the same path with an extra kick behind it.


#### But what about nodes?

Nodes are places with zero predicted particle hits.  In the double slit experiment, opening the second slit seems to place nodes where there were none with one split open.  In other words, locations that were perfectly valid with one slit open are suddenly avoided, even though that original slit is still unblocked.

But a node is infinitesimal in width.  Nodes or no, there are no locations of finite width with zero probability for a particle to hit.  What opening the second slit does, is allow for many new possible paths.  These new paths are not necessarily symmetrically distributed as far as where they land.  A somewhat likely detector location when one slit is open might not be so likely anymore if there are many new available paths and most of them don't go there.


#### But what about entanglement?

Entanglement is often portrayed in a one to one light, as if every particle that is a member of an entangled pair instantly knows about the other once measurement takes place.  However, to the best of my knowledge, the experimental and theoretical predictions are much more modest, merely a better than chance correlation.  

The conceptual model proposed here allows for enhanced correlations due to the induced clustering around integer units of action.  Two photons created from one event, for instance, are phase correlated.  They have their starting points at the same location, and on average both will move in 1 h units.  Thus, they maintain correlated probability distributions without needing to communicate.  The assumption is that other particle properties could also be jostled around, yet maintain the correlation as well.

Any particular measurement is chancy and uncertain, but the average over many measurement is well defined and correlated.


#### But what about the Uncertainty Principle?

In this model, the uncertainty principle is not so abstract, but a direct consequence of constant jostling.  A particle has an exact location and momentum at any time, but part of this momentum is due to the extra, unknown input from Brownian motion.  Several measurements must be taken to pin down the actual momentum.  But of course, this blows the position measurements.  


#### But what about the measurement paradox?

In this model, there is no measurement paradox, because the particle has a definite though chaotic history.  The wave function has no physical reality, and only represents a histogram.  Wave function collapse is no more physical then rolling a pair of dice while considering all possible outcomes.  Before the dice land, the histogram is the best basis for predicting the results.  But physically it is still particular dice in particular orientations the entire time.  The results appear to be guided by the histogram because after many rolls they will indeed do so.  But nobody would argue that before the dice landed they were just a histogram.  


#### But what about the observer effect?

Of course, the double slit experiment encompasses more than just interference patterns.  The interference goes away upon placing a detector at one slit.

But no detection occurs without interaction, likely with some exchange of momentum or subtle reflection.  Despite phase being a quantum property, the quantum level jostling makes no change to the phase on averaage.  It is only classical motion that actually changes the phase.  Think of the wave function for a bound state, such as the Hydrogen atom.  The phase only changes in the classical region around the nucleus where the kinetic energy is positive.  Far from the nucleus, the kinetic energy is negative, and the wave function is an ordinary decaying exponential.  

Therefore, a classical interaction at the slit can disrupt the phase relationships, smearing out the interference.  This is why there are weak measurements that disturb the interference pattern less and strong measurements that disturb the pattern more.  Upcoming simulations are planned to test this effect.


### References

- QED: The Strange Theory of Matter and Light.  Richard Feynman.
- Quantum Chemistry, Donald A McQuarrie.
- Introduction to Quantum Mechanics.  David J. Griffiths.


### Running the script

The script is presently fairly basic and in a state of flux. However feel free to try it out if you like.  A simulation may be run from the commandline using the following:

```buildoutcfg
$ python historians.py directory electrons
```

where "directory" is the name of a folder for storing the results, and "electrons" is the number of total successfully detected electrons.  Results are saved periodically, and the simulation may be stopped and restarted at will.  


#### Thank you!