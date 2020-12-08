#!/bin/bash

## Run simulations, fits, and chains for defined noise levels
#SBATCH --job-name=mcmcrun

## Allocation Definition
#SBATCH --account=astro
#SBATCH --partition=astro

## Nodes
#SBATCH --nodes=1

## Walltime (15 hrs)
#SBATCH --time=15:00:00

## Memory per node
#SBATCH --mem=128G

## Working directory
#SBATCH --chdir=/gscratch/astro/blindor/ttv_ss/

## Turn on e-mail notification
#SBATCH --mail-type=ALL
#SBATCH --mail-user=blindor@uw.edu

## Export variables to batch job session
#SBATCH --export=ALL

## Julia
export PATH=/gscratch/astro/blindor/juliastuff/julia:$PATH

julia full_run.jl test 30 &> results/test.out &
julia bgrun_10.jl try01 15 &> results/10try01.out &
julia bgrun_10.jl try02 30 &> results/10try02.out &
julia bgrun_10.jl try03 45 &> results/10try03.out &
julia bgrun_10.jl try04 60 &> results/10try04.out &
julia bgrun_10.jl try05 75 &> results/10try05.out &
julia bgrun_10.jl try06 90 &> results/10try06.out &
julia bgrun_10.jl try07 105 &> results/10try07.out &
julia bgrun_10.jl try08 120 &> results/10try08.out &
julia bgrun_10.jl try09 135 &> results/10try09.out &

wait
