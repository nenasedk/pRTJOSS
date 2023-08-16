---
title: 'Atmospheric Retrievals with petitRADTRANS'
tags:
  - Python
  - astronomy
  - exoplanets
  - atmospheres
authors:
  - name: Evert Nasedkin
    orcid: 0000-0000-0000-0000
    equal-contrib: true
    corresponding: true
    affiliation: 1 # (Multiple affiliations must be quoted)
  - name: Paul Mollière
    equal-contrib: true # (This is how you can denote equal contributions between multiple authors)
    affiliation: 1
  - name: Doriann Blain
    equal-contrib: true 
    affiliation: 1
  - name: Eleonara Alei
    equal-contrib: false 
    affiliation: 2
  - name: Karan Molaverdikhani
    equal-contrib: false
    affiliation: 3
  - name: Francois Rozet
    equal-contrib: false
    affiliation: 4
  - name: Aaron David Schneider
    equal-contrib: false
    affiliation: "5,6"
  - name: Tomas Stolker
    equal-contrib: false
    affiliation: 7
  - name: Nick Wogan
    equal-contrib: false 
    affiliation: 8
  - name: Mantas Zilinskas
    equal-contrib: false 
    affiliation: 7



affiliations:
 - name: Max Planck Institut für Astronomie, DE
   index: 1
 - name: ETH Zürich, Zürich, CH
   index: 2
 - name: Leiden Observatory, Leiden, NL
   index: 7
 - name: University of Washington, Seattle, USA
   index: 8
 - name: Ludwig-Maximilians-Universität, München, DE
   index: 3
 - name: University of Liège, Liège, BE
   index: 4
 - name: Centre for ExoLife Sciences, Niels Bohr Institute, Copenhagen, DK
   index: 5
 - name: Institute of Astronomy, KU Leuven, Leuven, BE
   index: 6
date: 14 August 2023
bibliography: paper.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
# aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
# aas-journal: Astrophysical Journal <- The name of the AAS journal.
#\footnote{\url{https://petitradtrans.readthedocs.io/en/latest/content/notebooks/pRT_Retrieval_Example.html}}. 
---

# Summary

`petitRADTRANS` (pRT) is a fast radiative transfer code used for computing emission and transmission spectra of exoplanet atmospheres [@molliere2019],
combining a FORTRAN back end with a Python based user interface.
These spectra can be used as a forward model for fitting spectroscopic data using Monte Carlo techniques, commonly referred to as an atmospheric retrieval [@madu2009].
Implementing such a framework in pRT is critical to remain at the forefront of exoplanet atmospheric science.
The retrieval module of pRT combines fast pRT with the `MultiNest` [@feroz2008; @feroz2009; @feroz2013; @buchner2014] and `Ultranest` [@buchner2014; @buchner2019] nested sampling codes, allowing for fast atmospheric retrievals on a large range of different types of exoplanet data.
Both samplers also offer MPI implementations, allowing for easy parallelisation across a large cluster.
With these new retrieval tools, it is now possible to use pRT to easily and quickly infer the atmospheric properties of exoplanets in both transmission and thermal emission.

# Statement of need

Atmospheric retrievals are a cornerstone of exoplanet atmospheric characterisation.
pRT is a powerful and user-friendly tool, and is unique in its abilities to characterise exoplanets in both emission and transmission.
Various thermal structures, chemistry and cloud parameterisations and opacity calculation methods can be combined and used to statistically identify the most favoured model for a given atmosphere.
With increasing volumes of both ground- and space-based spectra available, it is necessary for exoplanet researchers to have access to a range of characterisation tools.

# petitRADTANS Retrieval Module
The `Retrieval` module combines the `RadTrans` forward modelling class with a nested sampler via a likelihood function to perform an atmospheric retrieval.
Datasets, priors and other retrieval hyper parameters are set through the `RetrievalConfig` class, while the `models` module includes a range of complete atmospheric models that can be fit to the data.
Users can also define their own model function, making use of temperature profiles from the `physics` module and chemistry parameterisations from the `chemistry` module. 

Multiple datasets can be included into a single retrieval, with each dataset receiving its own `RadTrans` object used for the radiative transfer calculation, allowing for highly flexible retrievals where multiple spectral resolutions, wavelength ranges and even atmospheric models can be combined in a single retrieval.
Each dataset can also receive scaling factors (for the flux, uncertainties or both), error inflation factors and offsets.
The model functions are used to compute a spectrum $\vec{S}$, which is convolved to the instrumental resolution and binned to the wavelength bins of the data using a custom binning function to account for non-uniform bin sizes.
The resulting spectrum compared to the data with flux $\vec{F}$ and covariance $\mathbf{C}$ in the likelihood function:
\begin{equation}\label{eqn:loglike}
    -2\log\mathcal{L} = \left(\vec{S}-\vec{F}\right)^{T}\mathbf{C}^{-1}\left(\vec{S}-\vec{F}\right) + \log\left(2\pi\det\left(\mathbf{C}\right)\right).
\end{equation}
The second term is included in the likelihood to allow for uncertainties to vary as a free parameter during the retrieval, and penalizes overly large uncertainties.

pRT can compute spectra either using line-by-line calculations, or using correlated-k (c-k) tables for defining the opacities of molecular species.
We include up-to-date correlated-k line lists from Exomol [@tennyson2012; @mckemmish2016; @polyansky2018; @chubb2020] and HITEMP [@rothman2010], with the full set of available opacities listed in the online documentation.
The \verb|exo-k| package is used to resample the the correlated-k opacity tables to a lower spectral resolution in order to reduce the computation time [@leconte2021].
Combining the c-k opacities of multiple species requires mixing the distributions in $g$ space. 
Previously, this was accomplished by taking 1000 samples of each distribution.
This sampling process resulted in non-deterministic spectral calculations, leading to unexpected behaviour from the nested sampling process as the same set of parameters could result in different log-likelihood.
This has been updated to fully mix the c-k distributions.
Considering the first species, the second species is added in, and the resulting grid is sorted. 
The cumulative opacity grid is then mixed with the next species, a process which iterates until every species with significant opacity contributions (>0.1$\%$ of the current opacity in any bin) is mixed in to the opacity grid. 
Once complete, the resulting grid is linearly interpolated back to the 16 $g$ points at each pressure and frequency bin as required by pRT.
This fully deterministic process stabilized the log-likelihood calculations in the retrievals, and resulted in a 5$\times$ improvement in the speed of the c-k mixing function.

Various thermal, chemical and cloud parameterisations are available in pRT.
Built in temperature profiles range from interpolated splines to physically motivated profiles as in @guillot2010 and @molliere2020.
Equilibrium and disequilibrium chemistry can be interpolated from a pre-computed grid on-the-fly.
Chemical abundances can also be freely retrieved, with the additional possibility of using a combination of free and chemically consistent abundances.
Cloud parametersiations range from a 'grey' continuum opacity applied at all wavelengths, to clouds parameterised as in @ackerman2001, using log-normal or @hansen1971 particle size distributions with real optical opacities for different compositions and particle shapes, and including self-scattering.
Included in pRT is an option to use an adaptive pressure grid with a higher resolution around the location of the cloud base, and a lower resolution elsewhere. 
The higher resolution grid is 10 times as fine as the remaining grid, and replaces one grid cell above and below the cloud base layer, as well as the cloud base layer cell itself. 
This allows for more precise positioning of the cloud layers within the atmosphere. 

Photometric data are fully incorporated into the retrieval process.
The spectral model is multiplied by a filter transmission profile from the SVO database using the `species` package [@stolker2020].
This results in accurate synthetic photometry, which can be compared to the values specied by the user with the `add_photometry` function.

Publication ready summary plots of best fits, temperature and abundance profiles and corner plots can be automatically generated.
Multiple retrieval results can be combined in the plots for model intercomparisons.
Such results have been benchmarked against other widely used retrieval codes, in particular as part of the JWST Early Release Science program (Wellbanks et al, in prep).


# Acknowledgements

We acknowledge contributions 

# References
