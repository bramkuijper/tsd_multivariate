# Seasonality TSD: timing
Individual-based models on the evolution of temperature-dependent sex
determination in a seasonal environment in which seasonality coevolves
with sex allocation traits. This is all done in a patch-structured population
(linked by juvenile dispersal) with overlapping generations. 

## The environment
Temperature in each patch varies seasonally according to an intercept +
sinusoidal function + a patch-specific error sampled from a normal
distribuiton. The sinusoidal is given by
`sin(2 * pi * t / max_t)` where `t` is the current timestep of the simulation and
`max_t` is the duration of a single season. See `TSDSeasonal:update_environment()`
for more detail.

## Sex-specific consequences of the environment
We assume that optimal temperatures for hatchling survival are sex-specific,
where we have stabilizing selection around an optimum given by `theta[female]`
and `theta[male]`. For example, we could have `theta[female]=-0.5` and `theta[male]=0.5`,
and we can use inverse trig to calculate for which values of t those timepoints apply.

## Evolving loci
 * `t` (day of season): males and females are only able to reproduce
whenever `t` matches the current time step. 
This is governed by an evolving locus `t` that can vary from 0 to the length of the season.
 * `a`: reflects the pivotal
temperature (the temperature at which sex determination changes from producing a majority of 
one sex to producing a majority of the opposite sex)
 * `b` reflects the sensitivity
to temperature in sex determination. If `b=0`, then we have genetic sex determination, if
`b<0` then more females will be produced in colder envts, and more males will be produced
in warmer envts. If `b>0` then more males will be
produced in colder envts than females and more females will be produced in warmer envts.

## The question
Do individuals start to change the time of breeding if evnt is getting warmer. How does
TSD coevolve with that?



