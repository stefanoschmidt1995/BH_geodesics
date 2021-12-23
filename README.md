# BH motion

This simple script computes the motion of a satellite (test particle) around a BH. The computation can be performed both in GR or within the Netwonian limit.

## How it works

The geodesics around a non-spinning BH (Schwarzschild) are a simple matter: it's just about writing down the conservation of energy and integrating the trajectories. The [wikipedia page](https://en.wikipedia.org/wiki/Schwarzschild_geodesics) will tell you a lot on this.

## How to use it

The code works in python. You will need to be able to run python (better from a terminal) with `numpy`, `scipy` and `matplotlib` installed.

To compute the orbit, you need to write some options in a config_file.ini. Then you just type:

```Bash

    python BH_geodesics.py config_file.ini

```

A standard ini file looks like:
```

[BBH]
L= 5.02
phi_0= 0.
r_0 = 50.9
r_dot_0 = 0.
GR = 1
t_max= 9e4

max_step = 1e3
relative_units = 1
compare = 0
show = 1
save_trajectories = 0

```

Everything is in units of mass (where G = c = 1)
In the ini file, you need to set:

* `L` : angular momentum (in units of M**2)
* `r_0` : initial distance (in units of M)
* `phi_0` : initial angle (dimension less)
* `r_dot_0` : initial radial velocity (dimensionless)
* `GR`: a boolean variable to control whether the Newtonian potential or the GR potential shall be used
* `t_max` : maximum integration time (in units of M)

There are also a bunch of optionals parameters:

* `max_step` (default = 1e4): controls the maximum integration step: if the integration doesn't converge, it may be a good idea to decrase it 
* `relatives_units` (default = 0): if True, `r_dot_0` will be in units of sqrt(2*V(r,L)). This makes the dynamics more redable:
	- if `|r_dot_0|>1` there are unbound orbits (hyperbolic)
	- if `|r_dot_0|<1` we have a bound system (ellyptic)
* `compare` (default = 0): if set to True, both GR and newtonian solutions will be computed and plotted together
* `show`(default = 1): if True, it will plot show some plots
* `save_trajectories` (default = 0): if True, it will save the trajectories to a jpeg file

If you want to display a quick help message, just type `python BH_geodesics.py --help`.

You can also set several configurations in the same file. Make sure each sets of configurations starts with `[SET_NAME]`.

## Examples

In folder `plots`, you can find an ini file with some examples. Feel free to play with them and explore the different possibilities.

Setting a zero initial velocity and negative energy, you will find strongly precessing orbits: the trajectory will look very nice.

<p align="center">
  <img src="https://github.com/stefanoschmidt1995/BH_geodesics/blob/master/plots/nice_elliptic_trajectory.jpeg">
</p>

If we set a non zero outward radial velocity (still keeping negative total energy) we will obtain strongly precessing, strongly eccentric orbits.

<p align="center">
  <img src="https://github.com/stefanoschmidt1995/BH_geodesics/blob/master/plots/elliptic_trajectory.jpeg">
</p>

We can also set a positive total energy to obtain an hyperbolic trajectory. Starting far away from the BH and setting a large angular momentum, the path of the object will be slightly deviated.

<p align="center">
  <img src="https://github.com/stefanoschmidt1995/BH_geodesics/blob/master/plots/hyperbolic_trajectory.jpeg">
</p>

If we start closer to the black hole, the gravitational deviation effect will be much much higher.

<p align="center">
  <img src="https://github.com/stefanoschmidt1995/BH_geodesics/blob/master/plots/strong_hyperbolic_trajectory.jpeg">
</p>

Note that in every case, the effect of Newtonian physics and those of GR are very different! This is only because a 1/r**3 term in the potential: cool, isn't it?

Of course, GR is not always important. If we take the case of the earth orbiting around the sun, the Newtonian physics works just fine:

<p align="center">
  <img src="https://github.com/stefanoschmidt1995/BH_geodesics/blob/master/plots/earth_sun_trajectory.jpeg">
</p>

## Author

This was part of the tutorial for a GR course at Utrecht Univerisity. If you want more information, please feel free to send me an email: [stefanoschmidt1995@gmail.com](mailto:stefanoschmidt1995@gmail.com)












