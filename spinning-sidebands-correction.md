# Spinning Sidebands Correction
[Back to the Main page](README.md)

It is essential for the package to work that you know what are Spinning side bands, and how are corrected, as this explains why it is important to know the magnetic field stregth of the NMR and the spinning rate applied (which is now information requested by the NMR functions).

When doing solid state 13 CP MAS NMR if the rate of sample spinning is less than the frequency of the chemical anisotropy, spinning side bands (SSB) may appear in CPMAS spectra in addition to the signal centered at the isotropic chemical shift. SSB are positioned at chemical shifts equal to the spinning frequency at both the right and left side of the centerband, regardless of the strength of the applied magnetic field. SSB are relatively intense for natural biopolymers and soil organic matter where the chemical structures are inherently anisotropic. High rates of magic angle spinning can reduce the sideband intensity, but it is typically not feasible to remove the sidebands from anisotropic solids.  

The effect of spinning rate on the number and intensities of SSB is shown in Fig. 1. The SSB can be easily recognized because their intensities decrease and their positions change as the spinning rate varies, whereas the position and intensity of the centerband remain unchanged. The SSB can be easily recognized because their intensities decrease and their positions change as the spinning rate varies, whereas the position and
intensity of the centerband remain unchanged.

As shown below:

![image](https://github.com/LuisCol8/SOMnmR/assets/35764330/fd38b431-e32d-4b7e-be2f-5fbe2a7fb5ec)

Fig 1. Effect of fast rotation on spinning side bands (SSB). Increasing the sample spinning rate leads to decreasing of both amount and intensity of SSB. Ref: Duer M. J. (2002)

Different approaches have been suggested to suppress side bands. In particular, one of the best ways to eliminate SSB is to maximize the rotor spin velocity, i.e. up to 15 – 20 kHz.

Nevertheless, when available NMR devices do not allow fast spinning, other options may be applied. For instance, the method of Total Suppression of Side Bands (TOSS) pulse sequence. The TOSS sequence consists of a series of four to six 1808 pulses, with phase cycling to compensate for pulse imperfections, applied on the carbon-13 just before spectrum acquisition.

When the TOSS sequence is not available because of instrument limitation, a simple mathematical SSB subtraction has been suggested to account for the anisotropic signal dispersion in spectra of NOM. SSBs in CPMAS 13C-NMR spectra of NOM are mainly represented by CSA of carboxyl groups, while the SSB arising from aromatic systems are in most cases low and negligible. Aliphatic chains in amorphous materials do not produce side bands. Signal areas can be corrected for these SSBs by measuring the area under each side band in the region where it is clearly visible, and subtracting this value from the area of the region where the side band is covered by other signals. Since different carbon groups contribute to SSBs, their measured areas must be also added to the centerband area to account for the transfer of signal to the SSBs. A limitation of such mathematical SSB subtraction may be due to an underestimation of hidden SSBs (such as those possibly given by aromatic moieties), and to the assumption that visible SSBs have the same shape and intensity as the hidden ones.

Example of Mathematical substraction of SSB

If we asume that a magnetic field strength of 200 MHz, and at the spinning rate of 6.8 kHz is used, the Carbon-13 nucleus has a gyromagnetic ratio of ¼ relative to the 1H nucleus. 
Therefore, the 13C resonance frequency is 50 MHz (¼ * 300 MHz). A prerequisite for quantitative interpretation of NMR signal intensities (peak areas) is a magic angle spinning rate that exceeds the 13C frequency range of 200 ppm. At a frequency of 50 MHz, this corresponds to a minimum magic angle spinning (MAS) rate of 5,000 Hz (½ * 200 ppm x 50 MHz).

Therefore, the MAS rate of 6.8 kHz is sufficient to move the spinning sidebands outside of the 13C spectrum, allowing for quantitative interpretation of the signal intensities. However, 6.8 kHz MAS is not sufficient to eliminate sidebands entirely. The sidebands simply appear at a distance of 136 ppm from the isotropic signal (6,800 Hz / 50,000,000 Hz = 136 x 10-6 Hz). This separation of the sidebands from the isotropic peak allowed us to integrate the sideband peak areas and numerically add the sideband intensity to the isotropic signal intensity, in a procedure commonly called “sideband correction”.

A visual example of an NMR integral table, divided into the regions and their respective sidebands is in Fig. 2, here we can see in blue the two rows that correspond to the main integral of the Alkyl region, and in red the two regions that correspond to the SSB below (-136 to -91 ppm). Observe that this region below does not overalp with other region, while the region above the main peak (136 to 181 ppm) overalps with the Aryl and the Carboxyl region. Thus to obtain the "real" main peak area, we have to take this SSB which does not overlap and multiply it by 2 (to account for the one that we can not correct for). We cal aso observe that at the main peak region (0 to 45 ppm) there is an overlap from the SSB of Aryl and Carboxyl, which should be removed. we do this by substracting the values of these SSB where they are not overlaping with anything (272 to 317 ppm). We do this for each region.

![image](https://github.com/LuisCol8/SOMnmR/assets/35764330/84ad9b9e-480c-49ed-ac17-ac74a4bb35dd)


[Back to the Main page](README.md)
