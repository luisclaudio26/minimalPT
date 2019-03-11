# My road to bidirectional pathtracing
## Part I. The Tao of Light Transport Simulation

_What you should learn and then forget_

### First words ###

  This is the first part of a series where I'll be sharing with you my work, problems,
  and solutions towards while implementing a basic bidirectional pathtracing.
  I started with a blank screen and ended up with a "low fat" rendering system with
  as few features as I could get, just like _Raytracing in One Weekend_: spheres-only,
  simple BRDFs, no acceleration structure, no textures, but which *works*. This was
  my decision since the beginning, as I noticed in previous experiences implementing
  raytracers with a flexible architecture (like PBRT) and lots of features eventually
  gets you stuck in the future because of bad decisions which are now spread over a
  giant codebase. In other words: if you (like me) are inexperienced with rendering
  systems, trying to do something awesome as a first attempt will eventually lead
  to failure. PBRT/Mitsuba are nice because they were designed by experienced people
  who starts coding with a good mental map of what a rendering system should look like.

  I decided to implement bidirectional pathtracing because I wanted to go a bit
  beyond the basic pathtracers we see in computer graphics classes. At the same time,
  my experience with Physically-based Rendering (www.pbr-book.org) showed me that
  I *personally* cannot learn while simply following the text and trying to reproduce
  code, mainly because PBRTv3 is a big, well engineered system and this makes code
  more complicated to understand for a beginner - for example, the basic pathtracing
  code is full of code to handle subsurface scattering/volumetric effects, next
  event estimation and so on. Thus, my (ambitious) decision was to implement a BDPT
  system without looking at reference implementations, nor tutorials (which do not
  exist on the internet) nor textbooks; I wanted to "rediscover" the problems, come
  up with solutions to them and only look at textbooks to discover the efficient/elegant
  way of solving them or if I got too stuck and could not advance (which fortunately
  did not happen). Obviously I gave myself the right to lookup formulas for things
  multiple importance sampling which are mathematically more sophisticated.

  The aftermath of this experience is that now I'm way more confident about the
  concepts and basis knowledge on LTS: I now know where to look for missing cosine
  and "over $\pi$" terms, I know when I need to convert pdfs from solid angle to
  surface area domain and how to do it, I am more comfortable reasoning about
  importance sampling, etc. Reading the PBR book became infinitely easier and
  more informative, as they present nice solutions to problems I already faced.
  Reading complex articles on light transport also became easier, as my foundational
  knowledge on LTS is more solid.

### Basic vocabulary: radiometric quantities

  Before starting my journey, I noticed that I was not really familiar with the
  radiometric quantities and concepts involved in light transport simulation, so
  before typing code I read the first sections of chapter 3 of Eric Veach's thesis
  where he gives a good explanation on this, which I'll try to summarize. LTS is,
  before

  Light is nothing but a bunch of photons bouncing around, eventually reaching
  our eyes. What we perceive as bright/dim is simply more/less photons reaching
  our eyes, while color is the manifestation of photons of different wavelenghts
  arriving. Not surprisingly, the ultimate goal of light transport simulation is
  to count how many photons arrive at a certain point, so we can convert this
  to a pixel color somehow. The problem here is that counting photons itself is
  impractical, given that zillions of them may arrive in a point at a given time
  lapse; thus, we talk about _energy_, measured in Joules, which is an indirect
  way of talking about number of photons: it is easier to say a purely red lamp
  emits 10 J of enery in one second than explicitly saying that 123812938 photons
  [calculate this exactly] of wavelength x were emitted. Also, given that we talk
  about so many photons of tiny energy, it is natural to treat energy as a
  continuous quantity rather than discrete [explain this further].

  Once we defined energy, we can talk about _power_, which is lazily defined as
  "energy over time". Power is used to talk about how much energy per second arrives
  at a surface, or leaves a light source, etc. In LTS power is the main quantity
  we're interested in: given that we deal with systems in *equilibrium* [definition],
  it makes little sense to talk about energy (unless we're measuring the total
  energy arriving at a given sensor opened for a time period - which gives us
  motion blur), but even so we need not to reason that much about energy, and we
  end simply sampling power in different time steps before integrating.

  Once we can talk about how energy flows along time, we can "factor out" another
  dimension and talk about how much power arrives at a given _area_. One can easily
  see that the total power ariving at a surface area of 1m² does not tell that much
  about how power is _distributed_ over this surface: is all the power ariving at
  a single portion of this surface? or is it evenly spread? the quantity we use
  to talk about power arriving at a particular point is _irradiance_. The following
  picture is very common in graphics: [show hemisphere]. At a given point, photons
  may arrive from any direction; irradiance measures how many of that photons arive
  at this single point, per second. Just like summing all "the powers" give us
  the total energy, "summing all the irradiances" (i.e. integrating) of all points
  of a continuous surface gives us the total power arriving at this domain.

  Finally, recall the figure above. Considering a single point, it is natural then
  to talk about the photons arriving in a single direction: this is what _radiance_
  means, and is by far the quantity we'll deal with the most. Summing the energy
  of all directions in a single point gives us the irradiance at that point; summing
  all the irradiances of a surface, gives us the total power arriving at it;
  summing the power at all time steps, gives us the total energy that arrived
  there. [mostrar quadrinho com a ordem das coisas]

  The mathematics for these quantities is more or less intuitive: recall again
  the figuve above. Let's say $L(\mathbf{p}, \omega)$ is the irradiance function
  that measures the power arriving at a point $\mathbf{p}$ and direction $\omega$ -
  in other words, it (indirectly) tells us how many photons per second arrives at
  $\mathbf{p}$ from direction $\omega$. If we wish to know the total power arriving
  $\mathbf{p}$, we need to consider _all_ the directions that power may arrive from
  and sum all of them; as we are in a continuous setup, "summing" things is
  actually integrating $L(\mathbf{p}, \omega)$ over all possible values of $\omega$ -
  which is the set of all directions in the hemisphere above $\mathbf{p}$ - CG
  people usually call it $\mathcal{H}^+$:

  \[ L(\mathbf{p}) = \int_{\mathcal{H}^+} L(\mathbf{p}, \omega) \; d\omega\]

  Pun intended: we integrated (condensed) all the "angular information" of the
  irradiance function in a given point - we cannot distinguish from where direction
  did power came from, but we know that the total power arriving at $\mathbf{p}$
  is $L(\mathbf{p})$ - the _irradiance_ function. We "overloaded" the function
  $L$ because it still talks about power, but one version talks about power at
  a given a point, the other about power at a given point and direction. Analogously,
  if one integrates the irradiance $L(\mathbf{p})$ over all points $\mathbf{p}$ of
  a surface, we get the total power arriving at the surface.

  On the units: power is measured in Joules per sec (a.k.a. Watts); irradiance
  if (not surprisingly) measured in W/m² (or any other area unit). Radiance is
  measured in (W/m²)/sr = $\frac{W}{m² \cdot sr}$ (steradian, which is a measure of solid angle). Think of
  solid angle being to regular, planar angles what area is to length; thus,
  steradians are to radians what m² is to meter. Check the wikipedia article for
  a more precise definition, but you'll see it's a way of talking about (continuous)
  sets of directions over the hemisphere of directions - perfect, then, to talk
  about radiance.

### The physics of taking a digital picture

A digital camera is a _sensor_ with lots of fancy stuff like lenses and diaphragms
around it. This sensor is essentially a device with many tiny cells - we call it
_sensor pixels_ - that somehow "count" how many photons arrive at it: the more
photons it counts, the brighter will be this pixel.

---
The sensor of modern cameras usually don't care about the wavelength of these
photons. The sensor itself, then, just sees in shades of gray, brighter and dimmer.
Colored pictures are taken by placing a _color filter array_ on top of the sensor
[^1], so we partition pixel sensors in groups of red, green and blue: we give in spatial
resolution in exchange for "wavelength" resolution. This is solved by interpolating
each of the channels (output by the camera as a .RAW file or something like this)
to get three full sized images.

[^1]: 3-CCD cameras do not this because they have one separate sensor for each color.

---

---
Analog cameras (and the human eye is included!) do also "count" photons in a
certain sense, but their output is not linear as a digital camera usually is;
i.e. doubling the number of photons in a digital camera will make the image twice
as brighter, but that's not the case for analog cameras. This is a HUGE
difference between what a digital camera sees and what the human eye sees: our
eyes have a much larger _dynamic range_, and thus digital cameras are more limited.
Display monitors also suffer from this limitation when showing images.

Lots of research has been made to cope with this limitation of digital photography.
We'll also need to deal with this because the output of our renderer is an image
with pixels outside the usual [0,255] range and we'll need to find a way of converting
this to a range which our displays can show.

---

When we do physically-based rendering, we want to simulate this mechanism in order to achieve
an image which would be the same if taken with a camera in the real-world: indeed,
many heavy-duty renderers from the industry do support lenses and simulations of
real cameras, so it is easier to compose real footage and VFX.

When you have a finite aperture and a lens, each point in each pixel of the sensor
receives a stream of photons from a set of directions defined by the lens, the
aperture and their relative position to the point (see figure below).

---
When you have a pinhole camera, each point receives a stream of photons from
exactly one direction; just imagine what happens to the cone of light above when
it becomes narrower and narrower.

---

The whole rendering thing, then, is: _how do I count the number of photons (measure
the total energy) arriving at each pixel of this virtual sensor?_

Recall all the integrals and the radiometric quantities of the last section: we
are going to loop over each pixel of the sensor, integrating the power arriving
at the pixel over time. We'll compute this power by integrating the irradiance over
the square domain of the pixel sensor; finally, in order to compute irradiance,
we'll integrate the radiance arriving at the point from the set of directions
defined by the lens:

\[ P_k = \int_{t_s}^{t_e} \int_{A_k} \int_{A_l} L(s, s \leftarrow l, t) \; dl \, ds \, dt\]

Where:

$t_s$ and $t_e$ are the start and end times that the camera shutter was opened
$P_k$ is the output value of the k-th pixel (what we'll write to the final image)
$A_k$ is the set of points in the area of the k-th pixel
$A_l$ is the set of points in the area of the lens
$s \rightarrow l$ is a notation for the direction (the unit vector) that goes from $s$ to $l$
$L(s, \omega, t)$ is the radiance arriving at $s$ from direction $\omega$ at time $t$

I will not talk about integrating over time (which is the correct thing to do)
because I have no support for animations which would lead to motion blur. The
mathematical way of talking about this is considering that the radiance function
$L(s, \omega, t)$ is a _delta_ in the temporal dimension:

\[L(s, \omega, t) = L(s, \omega, t_0) \delta(t_0) \]

When you plug this into the integral above, it becomes:

\[ P_k = \int_{t_s}^{t_e} \int_{A_k} \int_{A_l} L(s, s \leftarrow l, t_0) \delta(t_0) \; dl \, ds \, dt = \int_{A_k} \int_{A_l} L(s, s \leftarrow l, t_0) \; dl \, ds\]

As we are always talking about the radiance at the same time instant $t_0$, we can
just drop it from the function and make it $L(s, \omega)$. Deltas are weird things
which deserve their own post so things won't look magical.

To summarize: you don't even need to remember that there was an integral over
time in first place, but you should know why we can forget this.

This, however, does not gives us a complete algorithm. First of all, we do not
know how to compute $L(s, \omega)$ - we do not even know what it depends on -,
nor we know how to compute these integrals which are so specific, defined over
weird domains.

Let's defer the discussion on $L(s,\omega)$ a bit and talk about weird integrals.
Maybe the standard tools from calculus (Green's theorem, change of coordinate
systems, etc.) would allows us to analytically integrate things (even more when
CG people deal with triangles almost 100% of the time!), but this is impractical:
just imagine a scene like _San Miguel_ with lots of objects everywhere; a radiance
function cannot get more weird than that [^2], so we take a _numerical_ approach.

[^2]: We can always think of Weierstrass and Dirichlet functions, which are indeed
worse (but they're not useful for generating beautiful images).

In this scenario, numerical integration is simpler, more flexible and allows us
to approximate stuff, which gives room to be faster then calculating exact values.
More importantly, numerical schemes for integration do not need much more info
than simply being able to evaluate the function in a given point (unlike analytical
integration, which "deep" knowledge about the domain). Finally, among the world
of numerical integration schemas, we pick _Monte Carlo integral_ because
it deals well with high dimension integrals and discontinuities, which is our
case.

---
Soon we'll see that $L(s,\omega)$ hides a HUGE integral inside it.

---


### Monte Carlo and its integrals

First of all, Monte Carlo integration is a whole complex subject in itself. Thus,
I'll give you a readers digest for this [awesome excerpt](https://cs.dartmouth.edu/wjarosz/publications/dissertation/appendixA.pdf) from prof. Wojciech Jarosz's[^3]
dissertation, which is more than enough to implement a simple Monte Carlo renderer
like ours.

[^3]: As far as I know, this is a polish name and is pronounced like *Voy-chair
Ya-rosh*, where the R in _chair_ is the same as in german "Achtung". I think
english speakers approximate this as *Voy-tek Ya-ross*.

First of all, recall that when we're dealing with a phenomenom/machine/algorithm
that spits random values, the mathematical way of talking about this is by using
_random variables_. A random variable $X$, thus, is the mathematical way of saying
"something that may be any value in a known range", and then this allows us to
operate with this "something" - doing things like Y = 2X+3, which is another random
variable. We define some functions that help us understand $X$: its expected value
$E[X]$, its variance $V[X]$, its cumulative distribution function $cdf(X)$ and probability
density function $pdf(X)$, for example. Properly talking about these things is beyond
the scope of this post, but this is described in any fundamental book on probability
(and on the reference linked above).

Monte Carlo integral is a tool that allows us to approximate the integral of a
function $f : \mathbb{R}^k \rightarrow \mathbb{R}$ by taking random samples on its domain and
computing an _estimator_:

<!-- \[ F = \frac{1}{N} \sum_{i = 1}^{N} \frac{f(X_i)}{pdf(X_i)} \] -->
\[ F = \frac{f(X_i)}{pdf(X_i)} \]

where each $X_i$ is a _random variable_ that "picks" a point on the domain $\mathbb{R}^k$,
and $pdf(X_i)$ is the _probability density function_ of this point, which intuitively
measures how likely it is to pick such a point in relation to the rest of the domain.

$F$ being a function of many random variables, it is a random variable itself. The
trick here is that the _expected value_ of $F$ is the integral of $f(x)$:

\[ E[ \langle F \rangle ] = E[ \frac{f(X_i)}{pdf(X_i)} ] = \int_{a}^{b} \frac{f(x)}{pdf(x)} pdf(x) \; dx = \int_{a}^{b} f(x) \; dx \]

So what? Well, the _strong law of large numbers_ tells us that, if we average a
bunch of samples from $F$, we converge to its expected value. In an algorithmic
fashion: we generate a bunch of random numbers $\xi_i$ with their PDFs, and take
the average of all $\frac{f(\xi_i)}{pdf(\xi_i)}$, we get _arbitrarily_ close to the integral we are looking for.

Notice how this fits our problem: we have a function $L(s, \omega)$ which takes
values on a set of points $s$ on some rectangle (the pixel in the sensor) and unit
vectors $\omega$ in 3D space, returning a radiance value. [^4] We are able to
evaluate $L$ on any $s$ and $\omega$, and being $s$ and $\omega$ simply points on
space, it is not hard to sample them. We do this a bunch of times and we can compute
the total power arriving at one pixel!

Finally, notice that the combination of a point $s$ and a direction $\omega$ defines
a ray in space - thus the name raytracing.

[^4]: This means we can consider $L(s, \omega)$ to be a function mapping values of some
subset $X \subset \mathbb{R}^6$ to values of $\mathbb{R}^+$.


---
I personally find PDFs hard to understand. They are defined as the derivative of
the _cumulative distribution function_ $CDF(X = k)$ of a random variable $X$, which measures
the probability of $X$ picking a value less than or equal to $k$; CDFs are monotonic,
non-negative, and always less than one. PDF's, however, can assume any non-negative
value as long as they integrate to one - they can even be infinity at one point
and zero in the rest, again the infamous delta function.

The notion of "relative likelyhood", however, matches well the PDF. For example:
think of an uniformly distributed variable - where all points in the domain are
equally likely to be selected. The CDF of two distinct points will be different,
but the PDF of any point on the domain is the same - as the CDF in this case is
a linear function, its derivative is a constant, $\frac{1}{b-a}$, were $a$ and
$b$ are the bounds of the 1D domain of $X$. In the context of Monte Carlo integration,
the PDF gives a lower _importance_ for samples which are more likely to be picked.

---

---
I'm not sure about how Monte Carlo integral works on non-euclidean domains nor on vector-valued
functions, so I prefer not to speculate about general domains and codomains. I'll update this
note as soon as I find something on this.

---


### Computing radiance: the Light Transport Equation

blabla

### Forgetting things (the right way)

As I told you in the beginning, this text is rather repetitive for a good reason:
in order to properly grasp the concepts taken from radiometry, one needs to reason
about the same thing multiple times, from different perspectives. To actually
understand things, I needed to read the same definition over and over, think of
it as departing from the definition of radiance, then do the same way departing
from the concept of energy, each time parsing the words "irradiance" as "power
per area per solid angle", then parsing "power" as "energy per second", then
thinking of "energy" as counting photons and so on.

blabla Taoism
