sckh
====

Programs for calculation of XAS and XES including vibrational effects

The main development is done in the src directory, and the code pieces will be moved there evenutally. 
For the time being specialized working programs can be found in the util directory.

Working codes (in util directory):

KH    - quantum mechanical solution of the Kramers-Heisenberg equations for vibrational effects in XES
SCKH  - semiclassical solution of the KH equations for XES
*sinc_DVR  - solution of 1d vibrational problem on a grid with a sinc DVR
*DVR_PO  - potential optimized DV
*Kramers-Heisenberg_resonant - like KH but for the resonant case (in the SCKH_resonant folder)

* these codes lack a logical build sequence right now...

Codes under development (in util directory):

XAS_nonadiabatic  - XAS and XAS including nonadiabatic effects. Also attempt to unify codes into one
SCKH_resonant  - like SCKH but for the reoannt case






