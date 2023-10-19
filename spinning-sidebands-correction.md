# Spinning Sidebands Correction
[Back to the Main page](README.md)

When doing solid state 13 CP MAS NMR if the rate of sample spinning is less than the frequency of the chemical anisotropy, spinning side bands (SSB) may appear in CPMAS spectra in addition to the signal centered at the isotropic chemical shift. SSB are positioned at chemical shifts equal to the spinning frequency at both the right and left side of the centerband, regardless of the strength of the applied magnetic field. SSB are relatively intense for natural biopolymers and soil organic matter where the chemical structures are inherently anisotropic. High rates of magic angle spinning can reduce the sideband intensity, but it is typically not feasible to remove the sidebands from anisotropic solids.  

The effect of spinning rate on the number and intensities of SSB is shown in Fig. 5. The SSB can be easily recognized because their intensities decrease and their positions change as the spinning rate varies, whereas the position and intensity of the centerband remain unchanged. The SSB can be easily recognized because their intensities decrease and their positions change as the spinning rate varies, whereas the position and
intensity of the centerband remain unchanged.

As shown below:

![image](https://github.com/LuisCol8/SOMnmR/assets/35764330/fd38b431-e32d-4b7e-be2f-5fbe2a7fb5ec)

Fig 1.Effect of fast rotation on spinning side bands (SSB). Increasing the sample spinning rate leads to decreasing of both amount and intensity of SSB. Ref: Duer M. J. (2002)

Example

If we asume that a magnetic field strength of 200 MHz, and at the spinning rate of 6.8 kHz is used, the Carbon-13 nucleus has a gyromagnetic ratio of ¼ relative to the 1H nucleus. Therefore, the 13C resonance frequency is 50 MHz (¼ * 300 MHz). A prerequisite for quantitative interpretation of NMR signal intensities (peak areas) is a magic angle spinning rate that exceeds the 13C frequency range of 200 ppm. At a frequency of 50 MHz, this corresponds to a minimum magic angle spinning (MAS) rate of 5,000 Hz (½ * 200 ppm x 50 MHz).

Therefore, the MAS rate of 6.8 kHz is sufficient to move the spinning sidebands outside of the 13C spectrum, allowing for quantitative interpretation of the signal intensities. However, 6.8 kHz MAS is not sufficient to eliminate sidebands entirely. The sidebands simply appear at a distance of 136 ppm from the isotropic signal (6,800 Hz / 50,000,000 Hz = 136 x 10-6 Hz). This separation of the sidebands from the isotropic peak allowed us to integrate the sideband peak areas and numerically add the sideband intensity to the isotropic signal intensity, in a procedure commonly called “sideband correction”.


[Back to the Main page](README.md)
