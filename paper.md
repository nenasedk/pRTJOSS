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
  - name: Tomas Stolker
    equal-contrib: false
    affiliation: 3
  - name: Nick Wogan
    equal-contrib: false 
    affiliation: 4
  - name: Karan Molaverdikhani
    equal-contrib: false
    affiliation: 5
  - name: Mantas Zilinskas
    equal-contrib: false 
    affiliation: 6
  - name: Francois Rozet
    equal-contrib: false
    affiliation: 7

affiliations:
 - name: Max Planck Institut für Astronomie, Germany
   index: 1
date: 14 August 2023
bibliography: paper.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
# aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
# aas-journal: Astrophysical Journal <- The name of the AAS journal.
#\footnote{\url{https://petitradtrans.readthedocs.io/en/latest/content/notebooks/pRT_Retrieval_Example.html}}. 
---

# Summary
The pRT codebase has undergone significant updates since its initial publication in `@molliere2019` 
A retrieval module combining the pRT spectrum calculations with the ``MultiNest` `[@feroz2008; @feroz2009; feroz2013]`and `Ultranest` `[@buchner2014]` samplers has been included to streamline retrievals of exoplanet atmospheres in emission and transmission.



# Statement of need
Multiple datasets can be included into a single retrieval, with each dataset receiving its own `RadTrans` object used for the radiative transfer calculation, allowing for highly flexible retrievals where multiple spectral resolutions, wavelength ranges and even atmospheric models can be combined in a single retrieval.
Each dataset can also receive scaling factors (for the flux, uncertainties or both), error inflation factors and offsets.
Several atmospheric models are built into the `models` module, allowing for a wide range of P-T, cloud and chemistry parameterizations.
These models are used to compute a spectrum $\vec{S}$, which is convolved to the instrumental resolution and binned to the wavelength bins of the data using a custom binning function to account for non-uniform bin sizes.
The resulting spectrum compared to the data with flux $\vec{F}$ and covariance $\mathbf{C}$ in the likelihood function:
\begin{equation}\label{eqn:loglike}
    -2\log\mathcal{L} = \left(\vec{S}-\vec{F}\right)^{T}\mathbf{C}^{-1}\left(\vec{S}-\vec{F}\right) + \log\left(2\pi\det\left(\mathbf{C}\right)\right).
\end{equation}
The second term is included in the likelihood to allow for uncertainties to vary as a free parameter during the retrieval, and penalizes overly large uncertainties.

pRT can compute spectra either using line-by-line calculations, or using correlated-k tables for defining the opacities of molecular species.
We include up-to-date correlated-k line lists from Exomol `[@tennyson2012; mckemmish2016; polyansky2018; chubb2020]` and HITEMP `[@rothman2010]`, with the full set of available opacities listed in the online documentation.
The \verb|exo-k| package is used to resample the the correlated-k opacity tables to a lower spectral resolution in order to reduce the computation time `[@leconte2021]`.

Included in pRT is an option to use an adaptive pressure grid with a higher resolution around the location of the cloud base, and a lower resolution elsewhere. 
The higher resolution grid is 10 times as fine as the remaining grid, and replaces one grid cell above and below the cloud base layer, as well as the cloud base layer cell itself. 
This allows for more precise positioning of the cloud layers within the atmosphere. 
Including this adaptive mesh, our pressure grid contains a total of 154 layers when two cloud species are used, which is the standard grid used in this work. 

Finally, photometric data are fully incorporated into the retrieval process.
As with spectroscopic data, a model is computed using a user-defined function.
This model spectrum is then multiplied by a filter transmission profile from the SVO database using the `species` package `[@stolker2020]`.
This results in accurate synthetic photometry, which can be compared to the values specied by the user with the `add_photometry` function.

## Correlated-k Implementation
The correlated-k implementation was significantly improved in both accuracy and speed.
Combining the c-k opacities of multiple species requires mixing the distributions in $g$ space. 
{\color{red}Previously, this was accomplished by taking 1000 samples of each distribution.}
This sampling process resulted in non-deterministic spectral calculations, resulting in unexpected behaviour from the nested sampling process, as the same set of parameters could result in varying log-likelihood.
This has been updated to fully mix the c-k distributions.
Considering the first species, the second species is added in, and the resulting grid is sorted. 
The cumulative opacity grid is then mixed with the next species, a process which iterates until every species with significant opacity contributions (>0.1$\%$ of the current opacity in any bin) is mixed in to the opacity grid. 
Once complete, the resulting grid is linearly interpolated back to the 16 $g$ points at each pressure and frequency bin as required by pRT.
This fully deterministic process stabilized the log-likelihood calculations in the retrievals, and resulted in a 5$\times$ improvement in the speed of the c-k mixing function.

# Using the Hansen distribution with EDDYSED
The EddySED cloud model from `@ackermann2001` is ...

Typically, it  a log-normal particle size distribution is assumed where the geometric particle radius will vary throughout the atmosphere as a function of the vertical diffusion coefficient $K_{\rm ZZ}$ and the sedimentation fraction $f_{\rm SED}$.
Here, we will substitute the log-normal particle size distribution with the Hansen distribution, and will rederive the calculation for the particle radius as a function of $K_{\rm ZZ}$ and $f_{\rm SED}$.

We begin with a review of the EddySED model: the distribution of the number of particles as a function of particle radius, $n(r)$ is approximated as a log-normal distribution with width $\sigma_{g}$ and characteristic geometric radius $r_{g}$.
\begin{equation}
    n(r) = \frac{N}{r\sqrt{2\pi}\log\sigma_{g}}\exp\left(-\frac{\log^{2}\left(r/r_{g}\right)}{2\log^{2}\sigma_{g}}\right),
\end{equation}
$N$ is the total number of cloud particles.

The goal of the EddySED model is to calculate $r_{g}$ for each layer in the atmosphere, given $K_{\rm ZZ}$ and $f_{\rm SED}$. 
It balances the upwards vertical mixing, parameterised by $K_{\rm ZZ}$ and the particle settling velocity, $v_{f}$
\begin{equation}\label{eqn:vf}
    v_{f} = w_{*}\left(\frac{r}{r_{w}}\right)^{\alpha}.
\end{equation}
Here $w_{*}$ is the convective velocity scale. Note that $r_{w}\neq r_{g}$. $r_{w}$ is the radius at which the  particle settling velocity equals the convective velocity scale:
\begin{equation}
    w_{*} = \frac{K_{zz}}{L},
\end{equation}
where $L$ is the convective mixing length.
Since $w_{*}$ is known, and $v_{f}$ can be found analytically as in `@ackermann2001` and `@podolak2003`, and a linear fit can be used to find both $\alpha$ and $r_{w}$.

With both of these quantities known, we follow AM01 and define $f_{\rm SED}$ as:
\begin{equation}\label{eqn:fsed}
    f_{sed} = \frac{\int_{0}^{\infty}r^{3+\alpha}n(r)dr}{r_{w}^{\alpha}\int_{0}^{\infty}r^{3}n(r)dr}
\end{equation}
For the log-normal distribution, one finds:
\begin{equation}
    \int_{0}^{\infty}r^{\beta}n(r)dr = Nr_{g}^{\beta}\exp\left(\frac{1}{2}\beta^{2}\log^{2}\sigma_{g}\right)
\end{equation}
Which we can then use to solve for $r_{g}$:
\begin{equation}
    r_{g} = r_{w}f_{sed}^{1/\alpha}\exp\left(-\frac{\alpha + 6}{2}\log^{2}\sigma_{g}\right)
\end{equation}

In order to use the Hansen distribution, we must recalculate the total number of particles $N$, and integrate the distribution for $f_{\rm SED}$. 
We note here that the Hansen distribution is parameterised by the effective radius, $\bar{r}$, rather than the geometric mean radius. 
In this derivation we do not correct for this difference in definition, as both act as nuisance parameters in the context of an atmospheric retrieval.

We start by giving the Hansen distribution in full:
\begin{equation}
    n(r) = \frac{N \left(\bar{r}v_{e}\right)^{\left(2v_{e}-1\right)/v_{e}}}{\Gamma\left(\left(1-2v_{e}\right)/v_{e}\right)} r^{(1-3v_{e})/v_{e}}\exp\left(-\frac{r}{\bar{r}v_{e}}\right)
\end{equation}
In `hansen1971` the authors use the parameters $a$ and $b$ to denote the mean effective radius and effective variance, which we write as $\bar{r}$ and $v_{e}$ respectively.
These differ from the simple mean radius and variance by weighting them by the particle area, as the cloud particle scatters an amount of light proportional to its area. Thus:
\begin{equation}
    \bar{r} = \frac{\int_{0}^{\infty}r\pi r^{2}n(r)dr}{\int_{0}^{\infty}\pi r^{2}n(r)dr}
\end{equation}
and 
\begin{equation}
    v_{e} = \frac{\int_{0}^{\infty} \left(r-\bar{r}\right)^{2} r^{2}n(r)dr}{\bar{r}^{2}\int_{0}^{\infty}\pi r^{2}n(r)dr}
\end{equation}

As in EddySED, we will fit for the settling velocity, which will provide us with $\alpha$ and $r_{w}$, which we can use to find $f_{\rm SED}$, as in \ref{eqn:fsed}.
However, we must now integrate the Hansen distribution. We find that:
\begin{equation}\label{eqn:hansint}
    \int_{0}^{\infty}r^\beta n_{Hans}(r)dr = \frac{v_{e}^{\beta} \left(v_{e}\beta + 2v_{e} + 1\right) \left(\frac{1}{\bar{r}}\right)^{-\beta} \Gamma\left(\beta + 1 + \frac{1}{v_{e}}\right)}{\left(-v_{e} + v_{e}^{\beta + 3} + 1\right) \Gamma\left(1 + \frac{1}{v_{e}}\right)}
\end{equation}
While this is complicated, when we can nevertheless use Eqns. \ref{eqn:fsed} and \ref{eqn:hansint} to solve for $\bar{r}$:
\begin{equation}
    \bar{r} = \left(\frac{ f_{sed}r_{w}^{\alpha}v_{e}^{-\alpha} \left(v_{e}^{3+\alpha} - v_{e} + 1\right) \Gamma\left(1 + \frac{1}{v_{e}}\right)}{\left(v_{e}\alpha + 2v_{e} + 1\right) \Gamma\left(\alpha + 1 + \frac{1}{v_{e}}\right)}\right)^{\frac{1}{\alpha}}.
\end{equation}
Thus for a given $K_{\rm ZZ}$, $f_{\rm SED}$ and $v_{e}$, we can find the effective particle radius for every layer in the atmosphere.

However, in order to compute the cloud opacity, we still require the total particle count. 
For a volume mixing ratio of a given species, $\chi_{i}$, we can integrate $n(r)$ to find $N$:
\begin{equation}
 N = \frac{\chi_{i}}{\left(\bar{r}^{3}v_{e} -1\right)\left(2v_{e} -1\right)}
\end{equation}

# Acknowledgements

We acknowledge contributions 

# References