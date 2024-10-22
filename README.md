Using the Solar System as a proxy for an exoplanetary system, 
given transit timing variations of Venus and Earth, carry out
a fit with TTVFaster and Nbody Grad to find period & mass of Jupiter.
Learn Julia the Hard Way:
https://github.com/chrisvoncsefalvay/learn-julia-the-hard-way/tree/master/_chapters

Julia Data Module:
https://github.com/JuliaIO/JLD2.jl
JPL ephemerides: https://ssd.jpl.nasa.gov/?planet_eph_export

DE440 : Created June 2020; compared to DE430, added about 7 years of new data.
        Referred to the International Celestial Reference Frame version 3.0.
        Covers JED 2287184.5, (1549 DEC 31) to JED 2688976.5, (2650 JAN 25)

Note: positions are integrated in astronomical units (AU), fixed AU = 149597870.700 km
but with polynomials stored in units of kilometers. 
The integration time units are days of barycentric dynamical time (TDB)

CALCEPH: https://github.com/JuliaAstro/CALCEPH.jl
retrieves the position, velocity and acceleration of Earth (geocenter) relative
to the Earth-Moon system barycenter in kilometers, kilometers per second and
kilometers per second square at JD= 2451624.5 TDB timescale 
For best accuracy the first time argument should be the integer part 
and the delta the fractional part (step through day).

compute(eph,jd0,time,target,center)

Compute position and velocity of target with respect to center at epoch
jd0+time. This method does not support the NAIF numbering scheme.
To get the best precision for the interpolation, the time is split in two
floating-point numbers. The argument jd0 should be an integer and time should
be a fraction of the day. But you may call this function with time=0 and jd0,
the desired time, if you don't care about precision.

Note: For planets without moons, Mercury and Venus, the barycenter location 
coincides with the body center of mass. compute(...) doesn't accept all NaifIDs
NaifID: 
    10          'SUN'
    0           'SOLAR SYSTEM BARYCENTER'
    2           'VENUS BARYCENTER'
    3           'EARTH MOON BARYCENTER'
    399         'EARTH'

Observer is calculating barycentic julian date for theings outside SS, with light travel time correction.

1). Simulate transit times from JPLEphemeris. Add noise option to data. 

sim_obs_and_find_times(jd1::Float64,sigma::Real,nyear::Real,obs::String)
NB:
12/18/23 : Somehow, routine with TTVFaster skips the transit for Venus (at t=14513.936223576813).
1/09/24 : Added break to while loop. Now stopping before t0[end]. pre-2024 runs dont have this constraint

2a). Carry out a linear fit to the transit times. 

2b). Call ttv_nplanet.jl from wrapper then compute the chi-square 
of the fit. Carry out an initial fit of the 2 inner planets
Note: jmax is the truncation of infinite [?] series in Laplace coefficients. It's larger when planets are closer to each other.
ttv_nplanet(nplanet::Int64,jmax::Int64,ntrans::Vector{Int64},data::Vector{T}) where T<:Real
ttv_wrapper(tt0,nplanet::Int64,ntrans::Vector{Int64},params::Vector{T},jmax::Integer,EM::Bool) 

fit_planet2(jd1::Float64,sigma::Real,nyear::Real,tref::Real,tol::Real,options::Array{String},save_as_jld2::Bool=false)

2c). Then add in the third (and later, fourth) planet. Initialize a grid of periods & 
phases of the outer planet, compute the best-fit at each.
Optimize the fit to the two sets of transit times by varying all of the
other parameters. 
    - subtracted tref from times to increase speed/accuracy. 
    - added tolerance limit to optimization 
    - tref, jd1, and tol defined in full_run.jl 
Note: This assumes 2 transits for Jupiter/Mars-analogues because that is required for TTVFaster
Show that likelihood curve peaks at real period of Jupiter.
Note: the third planet can't be too close to the transiting planets. 
    - Currently doing 5-22 years, but ran rebound simulations to test minimum orbit of jupiter analogue. Conditioned on retrieved V+E parameters, Jupiter orbit can go down to 1337.0 days but this doesn't account for the observed TTVs. Would need NbodyGradient to find actual minimum.

fit_planet3(filename::String,jd1::Float64,sigma::Real,nyear::Real,tref::Real,tol::Real,p3in::Float64,p3out::Float64,np3::Int,nphase::Int,obs::String)
do wide run first, for nper=500 points in grid
then do nper=200 around peak

2d). Search for moon. Initialize a grid of moon phases & compute the best-fit at each.
Optimize the fit to the two sets of transit times by varying all of the
other parameters. 
fit_moon(jd1::Float64,sigma::Real,nyear::Real,tref::Real,tol::Real,dpin::Float64,dpout::Float64,ndp::Int,nplanets::Real)

3).  Taking the minimum chi-square, or the maximum log-Probability,
run a markov chain with all 3 planets:  
-Assuming the host star has the same mass as the Sun, which 3 planets did you detect?  
-What are their masses and eccentricities (as well as uncertainties on these
quantities)?

3a). Run markov chain with planets + moon:
Note: deltaphi dist. is multi-modal, so ensure that value is inside dp range.

 MCMC(foutput::String,param::Array{Float64,1},lprob_best::Float64,
    nsteps::Int64,nwalkers::Int64,nplanet::Int64,ntrans::Array{Int64,1},
    tt0::Array{Float64,1},tt::Array{Float64,1},sigtt::Array{Float64,1},
    use_sigsys::Bool,EM::Bool) 

4).  Run routine for range of noise levels (sigmas) and observation time spans (nyears). 
Which runs best detect perturbing objects? 
For different time spans, there are diff. posterior distributions, 
must change grid search to accommodate for width of dist.

Creating contours in corner plots: take all pixels in 2d histogram, 
sort them from smallest to largest value. Find where 68% of values is about a point (i.e. greather than )

If number of independent samples for each walker >= 100, get results. 
09/22/2023 updated 12/13/2023
12/18/23 : Somehow, routine with TTVFaster skips the transit for Venus (at t=14513.936223576813). 
1/09/24 : Added break to while loop. Now stopping before t0[end]. Pre-2024 runs dont have this constraint, so they have   include a transit time after the end of our observation time (for each planet). Only noticed while comparing with NbodyGradient; and I'm not redoing the runs - since I don't think it would matter.
%The values in parentheses are the actual transit counts, while the preceding values were including a transit time after the end of our observation' time.  

##########################	Current State	##########################
0). Updated TTVFaster to be compatible with Julia v1.3+
1). With transit times of Earth & Venus, can constrain both of their masses
jd1=2.4332825e6; tref=2430000; tol=1e-5
1a). Can infer existence of Jupiter first then Moon
1b). clear gaussians in likelihood profiles, agrees with posterior dist.
2). - For 30 sec noise, ~18 years is enough to constrain Jupiter period from E+V 
        (mass?){
            - For 30 sec noise, ~22 years is enough to constrain Mars period fr}}'
    - For 60secs and 90 secs, 30 years is enough
    limits/constraints make sense based on TTVs
3). Less time span or more noise overestimate (or underestimate?) Jupiter period
        could be different definitions of period (linear ephemeris fit vs average time to orbit)
4). Once t_maxsinphi and t_maxcosphi are approx 0, no constrain on deltaphi
    Correlation betweem deltaphi and t_maxsinphi --> posterior broader than likelihood
5). Four planet model preferred over model with moon and 3 planets. 
        can't tell difference b/w Hppmp and Hpppp
6). If wrong four planet params, wrong Jupiter params due to long term Mars TTV signal.

Q). Degenaracy b/w Jupiter and Moon? <--tail on Jupiter period with moon
        if you don't have enough time spans
Q). P-M_p degeneracy? <--tail on Jupiter distributions
        assymetry for ecosw and esinw
Q). Wrong signs for evectors? <-- not first time this has been found

##########################  Project Tasks   ##########################
0). Complete bibliography. [  ]
0a). Find relevant papers and add them to .bib file [  ]
0b). Read and summarize relevant papers [ x ]
1). Write up model desctription (as above) [ x ]
2). Write up methods description [ x ]
3). Write up analysis description [ x ]
4). Create tables for parameters. [ x ]
5) Make likelihood profiles continuous [ ]
6). Make plots of orbits with 1-sigma uncertainties overplotted with the correct orbits. [  ] 
     - how to do this with eccentricities and omega? Need pomega and Omega?
7). Figure out whether the Earth-Moon barycenter offset causes
bias in measurements and if so, why.
8). Schedule parallel fit and chain runs on hyak.mox using slurm scheduler [  ]
     - exits prematurely, times out, or returns error in regress.jl assertion
10). Condense results to 1 equation fit (ex. how much of X to get Y uncertainty). [  ]

############################ Future Tasks   ##########################
Q: What really limits timing precision of Earth & Venus
about the Sun? 
1). How to account for missed transits in real data? 
multiple coefficeints -> loop over coefficients?
2). Could use existing telescope precision info 
3). Figure out what the actual expected timing precision would be (limited by stellar noise -- related to Tyler's work). 
Q: The masses inferred with sufficient data are good, although
still a bit more discrepant than I would like. need to implement an N-body fit?
4a). Refine TTVFaster estimates from finding Jupiter by applying NbodyGradient.
(should get better parameters for the masses of Venus and Earth)
Heirarchy example for Solar System:
Sun Venus Earth Moon Jupiter Saturn ....
indices = [[-1, 1, 0, 0, 0, 0],  # SUN & VENUS orbit in a binary
           [ 0, 0,-1, 1, 0, 0],  # EARTH & MOON orbit in a binary 
           [-1,-1, 1, 1, 0, 0],  # SUN & VENUS orbit about them 
           [-1,-1,-1,-1, 1, 0],  # (optional) Jupiter orbits about them
           [-1,-1,-1,-1,-1, 1],  # (optional) etc...
           [ 1, 1, 1, 1, 1, 1]]  # center of mass of the system
4c). Compare TTVFaster and NBody Grad fits. 
5a). Figure out how to speed things up so I can do a global
search, and explore duration & error bar dependence. 
5b). Do inverse matrix fitting for linear parameters (Jupiter period & Moon deltaphi) to speed things up (might be more robust).
5c). Maybe make a struct to hold the pre-computed Laplace coefficents,
and pass this to routines, or create a closure for this.
5d). Analytical Jacobian (could speed up TTVFaster). How do finite differences compare to analytic derivatives?
    - computate Laplace coefficients, need derivatives of coeffs from table 1 of Agol&Deck 2015
    - write tests to show that analytic derivatives work.
6). From posteriors, show how well we can measure mean insolation. 
7). Show that model is correct either way (Moon first then Jupiter). 
    - unrealistic because it's more likely that the giant planet would be discovered first since it's easier
8). Make model of actual transit light curves (as opposed to just transit times).
    - Show how well constrained densities are (for Earth and Venus).
    - Show how well constrained densities are for Sun.

   
<!-- 
##########################  Completed Tasks  ##########################
9). See if we can detect Mars [ x ] or Saturn. [ x ]
8). See which scenario best fits simulated data [ x ]
7). Add M_p > 0 prior to MCMC [ x ]
6). Add in 4th planet. Fit for best params [ x ]
6a). Search for second peak in likelihood profile to fit. [ x ]
5). Analyze chain results: trace plots, uncertainties, etc.
5a). See how many observations would be needed (minimum number of years required). [ x ]
5b). See what the necessary precision would be (vary noise added to simulations). [ x ]
4). Show models are correct: derived Earth and Venus parameters.
4b). Make plots of histograms of parameter results from MCMC with correct values. [ x ]
4d). Make plots of logL for Jupiter period and Moon deltaphi with correct values at peak. [ x ] 
    (include histograms of posterior results) 
4e). Make plots of posterior results of model fit to simulated times. [ x ] 
3). Create slurm file to run multiple grids on hyak.mox. [ x ] 
2a). Add in the option for Moon. [ x ]
2b). Fit for Moon deltaphi. [ x ]
1). Makes plots of the contributions of individual bodies [ x ]
    (including the ones we are neglecting).   


#### Hyak  slurm example
julia full_run.jl grid 30.0 40 ppp &> results/run.out &
julia EMB_run.jl grid 10 40 ppppp &> results/p5test.out &
......
obs = "fromEMB" or "fromEV"
label = [ppp, ppmp, pppp, etc.]
runtype = [sim, grid, mcmc, wide]
runtype, label = ARGS[1], ARGS[4]
sigma, nyear = parse(Float64,ARGS[2]),parse(Float64,ARGS[3])
np3 = [fine=200, medium=100, coarse=50]     # test=10 
nphase = [fine=72, medium=36, coarse=18]    # test=10 
ndp = [fine=180, medium=72, coarse=36]      # test=10 
steps=[short=10000, med=50000, long=100000] # test=1000
sigmas = [10, 30, 45, 60, 75, 90, 105, 120] # which of these are realistic?
years = [10, 12, 15, 18, 20, 23, 25, 28, 30, 40] # how often to check results?
......
####
Plotting examples
plot_contrib(30,30,["fromEMB","p4","best_p4"])

#### More general example for github
The data array contains parameters that describe a multi-transiting planet system. In the case of 2 planets, there are 10 parameters. 
```julia
  # Set up data structure to hold planet-plane properties,passed to TTVFaster
  data=init_param
  julia> p1=TTVFaster.Planet_plane_hk(data[1],data[2],data[3],data[4],data[ 5]);
  julia> p2=TTVFaster.Planet_plane_hk(data[6],data[7],data[8],data[9],data[10]);
  # Compute expected transit times (if there were no perturbations): 
  time1 = collect(p1.trans0 .+ range(0,stop=nt1-1,length=nt1) .* p1.period);
  time2 = collect(p2.trans0 .+ range(0,stop=nt2-1,length=nt2) .* p2.period);
  # Initialize the computation of the Laplace coefficients:
  ttv1 = zeros(nt1);
  ttv2 = zeros(nt2);
  # Need first call to TTVFaster,without optimizing
  julia> dummy=TTVFaster.compute_ttv!(jmax,p1,p2,time1,time2,ttv1,ttv2)
```
