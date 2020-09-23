Using the Solar System as a proxy for an exoplanetary system, 
given transit timing variations of Venus and Earth, carry out
a fit with TTVFaster and Nbody Grad to find period & mass of Jupiter.

CALCEPH 
retrieves the position, velocity and acceleration of Earth (geocenter) relative
to the Earth-Moon system barycenter in kilometers, kilometers per second and
kilometers per second square at JD= 2451624.5 TDB timescale

For best accuracy the first time argument should be the integer part 
and the delta the fractional part (step through day)

Note: For planets without moons, Mercury and Venus, 
the barycenter location coincides with the body center of mass. 

NaifID: 
      2           'VENUS BARYCENTER'
      3           'EARTH MOON BARYCENTER'
      399         'EARTH'

1). Simulate transit times from JPLEphemeris. Add noise option to data.

sim_times has multiple methods with the following arguments
sim_times(jd1::Float64, jd2::Float64, jdsize::Int64, 
    addnoise::Bool=false, sigma::Float64=0.0, EMB::Bool=true, seed::Int=42)

2a). Carry out a linear fit to the transit times. 

2b). Call ttv_nplanet.jl from wrapper then compute the chi-square 
of the fit. Carry out an initial fit of the 2 inner planets

2c). Then add in the third planet. Initialize a grid of periods & 
phases of the outer planet, compute the best-fit at each.
Optimize the fit to the two sets of transit times by varying all of the
other parameters. 

Note: This assumes 2 transits for Jupiter because transits are required for TTVFaster

Show that likelihood curve peaks at period of Jupiter.
Note: the third planet can't be too close to the transiting planets.

fit_mysteryplanet3(filename::String, label::String,
  p3in::Float64=4000.0, p3out::Float64=4600.0, np3::Int=10, nphase::Int=10, 
  addnoise::Bool=false, sigma::Float64=0.0, EMB::Bool=true)

3).  Taking the minimum chi-square, run a markov chain with
all 3 planets.  Assuming the host star has the same mass
as the Sun, which 3 planets did you detect?  What are their
masses and eccentricities (as well as uncertainties on these
quantities)?

MCMC(param::Array{Float64, 1},label::String,
  nsteps::Int64,nwalkers::Int64,nplanet::Int64,ntrans::Array{Int64, 1},
  tt0::Array{Float64, 1},tt::Array{Float64, 1},sigtt::Array{Float64, 1},
  EMB::Bool=true,use_sigsys::Bool=true)  

4). Refine TTVFaster estimates from finding Jupiter by applying NbodyGradient.
(should get better parameters for the masses of Venus and Earth)

Heirarchy example for Solar System:

Sun Venus Earth Moon Jupiter Saturn ....
indices = [[-1, 1, 0, 0, 0, 0],  # SUN & VENUS orbit in a binary
           [ 0, 0,-1, 1, 0, 0],  # EARTH & MOON orbit in a binary 
           [-1,-1, 1, 1, 0, 0],  # SUN & VENUS orbit about them 
           [-1,-1,-1,-1, 1, 0],	 # (optional) Jupiter orbits about them
           [-1,-1,-1,-1,-1, 1],	 # (optional) etc...
           [ 1, 1, 1, 1, 1, 1]]  # center of mass of the system

Learn Julia the Hard Way:
https://github.com/chrisvoncsefalvay/learn-julia-the-hard-way/tree/master/_chapters

Julia Data Module:
https://github.com/JuliaIO/JLD2.jl

can save fit_p3 data:
julia> using JLD2
julia> @save "OUTPUTS/p3_fit.jld2" param_p3
julia> @save "evj_mcmc_01_chi.jld2" chi_mcmc
julia> @save "evj_mcmc_01_par.jld2" par_mcmc

The, I can restore these later:
julia> using JLD2
julia> @load "OUTPUTS/p3_fit_test.jld2"

julia ttv_like_planet_b.jl &> ttv_likelihood_planetb_3.0sig.txt &

TODO:
-define format for grids and MCMC runs 
run_types = [extrashort=(4230-4430), short=(2000-5000), medium=(700-10000), wide=(500-18000)]
grid_types = [extrafine=1000, fine=100, medium=10, coarse=2]
noise = [10.0, 15.0, 30.0, 45.0, 60.0, 120.0, 240.0]
label 	run_type	grid_type	noise	
test 	extrashort	coarse		30
try001	short	 	medium		30
try002	medium		fine		30
try003	medium		fine

8/4/2020
##########################	Current State	##########################
0). Updated TTVFaster to be compatible with Julia v1.1
1). With transit times of Earth & Venus, can infer both of
their masses, as well as existence of Jupiter first then Moon
2). Q: What really limits timing precision of Earth & Venus
about the Sun? (related to Tyler's work)
3). The masses inferred with sufficient data are good, although
still a bit more discrepant than I would like:  I need to
implement an N-body fit.


##########################	Next Tasks	##########################
1). Makes plots of the contributions of individual bodies (including the ones we are neglecting). [ x ]
2a). Add in the option for Moon. [ x ]
2b). Fit for Moon deltaphi. [ x ]
3a). Create slurm file to run multiple grids on hyak.mox [ ]
3b). See how many observations would be needed (minimum number of years required). [ ]
3c). See what the necessary precision would be (add noise to simulations). [ ]
4a). Show models are correct: derived Earth and Venus parameters.
4b). Make plots of histograms of parameter results from MCMC with correct values. [ ]
4c). Make plots of orbits with 1-sigma uncertainties overplotted with the correct orbits. [ ]
4d). Make plots of logL for Jupiter period and Moon deltaphi with correct values at peak. [ ]
5a). Using TTVFaster for first estimate, do NBody Gradient fit. [ ]
5b). Compare TTVFaster and NBody Grad fits. [ ]
5c). Make plots of posterior results of model fit to simulated times. [ ] 
6). From posteriors, show how well we can measure mean insolation (eccentricity of Earth's orbit). [ ]
7). Show that model is correct either way (Moon first then Jupiter). [ ]


##########################	Writing Tasks	##########################
1). Write up model desctription (as above). [ ]
2a). Find relevant papers and add them to .bib file [ x ]
2b). Read and summarize relevant papers [ ]


##########################	Optional Tasks	##########################
3). Figure out what the actual expected timing precision
would be (limited by stellar noise -- related to Tyler's work). 
3a). Could use existing telescope precision info
5). See if we can detect Mars [ ] or Saturn.
8). Figure out why the Earth-Moon barycenter offset causes
bias in measurements (see if this is the case).
9a). Figure out how to speed things up so I can do a global
search, and explore duration & error bar dependence. 
9b). Do inverse matrix fitting for linear parameters (Jupiter period & Moon deltaphi) to speed things up (might be more robust).
9c). Maybe make a type to hold the pre-computed Laplace coefficents,
and pass this to routines, or create a closure for this.
13a). Make model of actual transit light curves.
13b). Show how well constrained densities are (for Earth and Venus).
13c). Show how well constrained densities are for Sun.

