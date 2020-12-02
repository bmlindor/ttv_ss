#!/bin/bash

## Job Name
#SBATCH --job-name=mcmcrun

## Allocation Definition
#SBATCH --account=astro
#SBATCH --partition=astro

## Nodes
#SBATCH --nodes=1

## Walltime (3 hrs)
#SBATCH --time=15:00:00

## Memory per node
#SBATCH --mem=128G

## Working directory
#SBATCH --chdir=/gscratch/astro/blindor/research/ttv_ss

## Turn on e-mail notification
#SBATCH --mail-type=ALL
#SBATCH --mail-user=blindor@uw.edu

## Julia
export PATH="/gscratch/astro/blindor/juliastuff/julia:${PATH}"
export JULIA_PKGDIR="/usr/lusers/blindor/.julia"
##export JULIA_DEPOT_PATH="/gscratch/astro/blindor/julia-1.0.5/bin/"

julia full_run.jl &> try001/


