# FluctReconPy

Electron density fluctuation reconstruction in BES measurements by minimalization of undulation.

The project has been unarchived because reconstruction of the individual channel fluctuations are better than for the direct method. The estimation of the blob sizes and positions were based on an averaging scheme which averaged out the noise in both cases. However, the method does work properly and the result is more convincing than with the direct method.

The method is demonstrated on the KSTAR viewing geometry but it could be applied to any other beam measurement at the edge of the plasma where light fluctuations can be taken linearly dependent on the density fluctuations.

The primary language of this module is IDL (despite the name FluctReconPy). A python version is available as well but it fails the benchmark agains the IDL version due to different and unresolved matrix multiplications and noise simulation.
