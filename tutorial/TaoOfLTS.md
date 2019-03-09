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
_sensor pixels_ - that somehow "count" how many photons arrive at it.

### Computing radiance: the Light Transport Equation

blabla

### Monte Carlo and its integrals

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
