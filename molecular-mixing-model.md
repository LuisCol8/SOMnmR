# Molecular Mixing Model
[Back to the Main page](README.md)

The molecular mixing model is a method to predict the molecular composition of a sample. It calculates the mathematically optimal composition of a sample using model components together so that the N:C ratio and the signal intensity in some regions of the NMR spectrum for the model mixture equalled those for the sample.
This proces was accomplished by solving Equations 1 to 6 (Fig 1.), the solution of these equations deﬁned a proportion of each model component for that sample (a to f). The proportions of the various model components are then used to calculate a predicted spectral distribution for that sample by summing, for each chemical shift region, the product of (a) the proportion ofeach
model component and (b) the present C in that region for that model component.

![image](https://github.com/LuisCol8/SOMnmR/assets/35764330/a99bd9b7-8bf2-4f1c-b806-cc155a22fe75)

Equations (1)–(6), a, b, c, d, e and f equal the proportions of components A (carbohydrate), B (protein), C (lignin), D (aliphatic material), E (carbonyl) and F (char) in the model. Equation (1) ensures that the sum of all component proportions equals 1. In Equation (2), n equals the N:C ratio of the component (or sample) speciﬁed (e.g., nA equals the N:C ratio of component A, the carbohydrate component). Molar N:C ratios ofthe model components are shown in Table 1. In Equations (3)–(6), a, b, v and d equal the proportions ofcarbon in the speciﬁed components (or sample) that resonate in the 45 to 	10, 95–60, 210–165 and 145–110 ppm chemical shift regions, respectively (e.g., aA equals the proportion of total signal in the 45 to 10 ppm region in the carbohydrate component).

The original paper from Nelson et al. (1999), used an amalgam of pure compounds as references. The  lignin is  assumed to be equivalent to  the average of those acquired for spruce and red alder lignin  presented by Wilson (1987), and the lipid is calculated from structures proposed by Kolattukudy (1980).

### Table 1. Reference compounds for natural soil organic matter, as published by Nelson et al. (1999).
| Chemical shift region (ppm)   | Aminoacid | Hexose  | Lignin | Cutin  | Suberin | Char  | Chitin | Carbonyl  |
| ------------- | ------------- | ------------- | ------------- | ------------- | ------------- | ------------- | ------------- | ------------- |
| 210 to 165 (Amide/Carboxyl) | 26.4 | 0.0  | 4.6  | 6.6  | 4.3  | 5.6  | 12.5  | 100.0  |
| 165 to 145 (Phenolic)  | 2.5  | 0.0  | 19.5  | 0.7 | 6.2  | 16.1  | 0.0  | 0.0  |
| 145 to 110 (Aromatic)  | 7.5  | 1.0  | 30.5  | 3.6  | 19.1  | 73.9  | 0.0  | 0.0  |
| 110 to 95 (Di-O-Alkyl)  | 0.0  | 15.7  | 8.6  | 0.0  | 0.0  | 4.3  | 12.5  | 0.0  |
| 95 to 60 (O-Alkyl)  | 2.1  | 79.0  | 12.5  | 9.0  | 6.8  | 0.0  | 50.0  | 0.0  |
| 60 to 45 (N-Alkyl/Methoxyl)  | 4.3  | 21.9  | 13.8 | 2.9  | 0.0  | 0.0  | 12.5  | 0.0  |
| 45 to -10 (Alkyl)  | 39.6 | 0.0  | 10.5  | 75.6  | 60.7  | 0.0  | 12.5  | 0.0  |
| Molar N:C   | 0.32  | 0.0  | 0.0  | 0.0  | 0.0  | 0.0  | 0.125  | 0.0  |

A subsequent paper from Baldock et al (2004), uses a modified version of this, adjusted for terrestrial, and aquatic ecosystems. Here it is considered for cabohydrates a cellulose model, for protein based on extractable aminoacids from Friedel and Scheller (2002).

### Table 2. Reference compounds for natural terrestrial organic matter, as published by Baldock et al (2004).
| Chemical shift region (ppm)   | Carbohydrate | Protein  | Lignin | Lipid  | Carbonyl | Char  |
| ------------- | ------------- | ------------- | ------------- | ------------- | ------------- | ------------- |
| 210 to 165 (Amide/Carboxyl) | 0.0 | 28.3  | 4.6  | 6.6  | 100  | 3.9  |
| 165 to 145 (Phenolic)  | 0.0  | 1.3  | 19.5  | 0.7 | 0.0  | 15.2  |
| 145 to 110 (Aromatic)  | 0.0  | 8.9  | 30.5  | 3.6  | 0.0  | 72.0  |
| 110 to 95 (Di-O-Alkyl)  | 16.7  | 0.0  | 8.6  | 0.0  | 0.0  | 5.3  |
| 95 to 60 (O-Alkyl)  | 83.3  | 3.5  | 12.5  | 9.0  | 0.0  | 1.8  |
| 60 to 45 (N-Alkyl/Methoxyl)  | 0.0  | 22.6  | 13.8 | 4.5  | 0.0  | 1.7  |
| 45 to -10 (Alkyl)  | 0.0 | 35.4  | 10.5  | 75.6  | 0.0  | 0.0  |
| Molar N:C   | 0.0  | 0.275  | 0.0  | 0.0  | 0.0  | 0.0  |
| C   | 1.0 | 1.0  | 1.0 | 1.0  | 1.0 | 1.0  |
| N   | 0.0 | 0.27  | 0.0 | 0.0  | 0.0 | 0.0  |
| H   | 1.67 | 1.1  | 1.24 | 1.94  | 1.0 | 0.45  |
| O   | 0.83 | 0.16  | 0.43 | 0.24  | 2.0 | 0.41  |


### Table 3. Reference compounds for natural aquatic organic matter, as published by Baldock et al (2004).
| Chemical shift region (ppm)   | Carbohydrate | Protein  | Lignin | Lipid  | Carbonyl | Char  |
| ------------- | ------------- | ------------- | ------------- | ------------- | ------------- | ------------- |
| 210 to 165 (Amide/Carboxyl) | 0.0 | 26.4  | 4.6  | 5.6  | 100  | 3.9  |
| 165 to 145 (Phenolic)  | 0.0  | 2.5  | 19.5  | 0.0 | 0.0  | 15.2  |
| 145 to 110 (Aromatic)  | 0.0  | 7.5  | 30.5  | 0.0  | 0.0  | 72.0  |
| 110 to 95 (Di-O-Alkyl)  | 16.7  | 0.0  | 8.6  | 0.0  | 0.0  | 5.3  |
| 95 to 60 (O-Alkyl)  | 83.3  |  2.1  | 12.5  | 0.0  | 0.0  | 1.8  |
| 60 to 45 (N-Alkyl/Methoxyl)  | 0.0  | 4.3  | 13.8 | 0.5  | 0.0  | 1.7  |
| 45 to -10 (Alkyl)  | 0.0 | 39.6   | 10.5  | 94.4  | 0.0  | 0.0  |
| Molar N:C   | 0.0  | 0.32 | 0.0  | 0.0  | 0.0  | 0.0  |
| C   | 1.0 | 1.0  | 1.0 | 1.0  | 1.0 | 1.0  |
| N   | 0.0 | 0.27  | 0.0 | 0.0  | 0.0 | 0.0  |
| H   | 1.67 | 1.17  | 1.24 | 1.89  | 1.0 | 0.45  |
| O   | 0.83 | 0.12  | 0.43 | 0.11  | 2.0 | 0.41  |

Later Nelson et al. (2005), published another set of reference compounds, which consider the carbohydrate and the protein fraction are based on the hexose and aminoacid composition described above and in Nelson (1999).

### Table 4. Reference compounds for natural terrestrial organic matter, as published by Nelson et al. (2005).
| Chemical shift region (ppm)   | Carbohydrate | Protein  | Lignin | Lipid  | Carbonyl | Char  |
| ------------- | ------------- | ------------- | ------------- | ------------- | ------------- | ------------- |
| 210 to 165 (Amide/Carboxyl) | 0.0 | 30.4  | 4.6  | 6.6  | 100  | 5.7  |
| 165 to 145 (Phenolic)  | 0.0  | 1.0  | 19.5  | 0.7 | 0.0  | 16.1  |
| 145 to 110 (Aromatic)  | 1.0  | 4.5  | 30.5  | 3.6  | 0.0  | 73.9  |
| 110 to 95 (Di-O-Alkyl)  | 15.7  | 0.0  | 8.6  | 0.0  | 0.0  | 4.3  |
| 95 to 60 (O-Alkyl)  | 79.0  | 2.9  | 12.5  | 9.0  | 0.0  | 0.0  |
| 60 to 45 (N-Alkyl/Methoxyl)  | 4.3  | 24.7  | 13.8 | 4.5  | 0.0  | 0.0  |
| 45 to -10 (Alkyl)  | 0.0 | 36.6  | 10.5  | 75.6  | 0.0  | 0.0  |
| Molar N:C   | 0.0  | 0.266  | 0.0  | 0.0  | 0.0  | 0.0  |
| C   | 1.0 | 1.0  | 1.0 | 1.0  | 1.0 | 1.0  |
| N   | 0.0 | 0.27  | 0.0 | 0.0  | 0.0 | 0.0  |
| H   | 1.67 | 1.1  | 1.24 | 1.94  | 1.0 | 0.45  |
| O   | 0.83 | 0.16  | 0.43 | 0.24  | 2.0 | 0.41  |


### Table 5. Reference compounds for aquatic terrestrial organic matter, as published by Nelson et al. (2005).
| Chemical shift region (ppm)   | Carbohydrate | Protein  | Lignin | Lipid  | Carbonyl | Char  |
| ------------- | ------------- | ------------- | ------------- | ------------- | ------------- | ------------- |
| 210 to 165 (Amide/Carboxyl) | 0.0 | 30.4  | 4.6  | 5.6  | 100  | 5.7  |
| 165 to 145 (Phenolic)  | 0.0  | 1.0  | 19.5  | 0.0 | 0.0  | 16.1  |
| 145 to 110 (Aromatic)  | 1.0  | 4.5  | 30.5  | 0.0  | 0.0  | 73.9  |
| 110 to 95 (Di-O-Alkyl)  | 15.7  | 0.0  | 8.6  | 0.0  | 0.0  | 4.3  |
| 95 to 60 (O-Alkyl)  | 79.0  | 2.9  | 12.5  | 0.0  | 0.0  | 0.0  |
| 60 to 45 (N-Alkyl/Methoxyl)  | 4.3  | 24.7  | 13.8 | 0.5  | 0.0  | 0.0  |
| 45 to -10 (Alkyl)  | 0.0 | 36.6  | 10.5  | 94.4  | 0.0  | 0.0  |
| Molar N:C   | 0.0  | 0.266  | 0.0  | 0.0  | 0.0  | 0.0  |
| C   | 1.0 | 1.0  | 1.0 | 1.0  | 1.0 | 1.0  |
| N   | 0.0 | 0.27  | 0.0 | 0.0  | 0.0 | 0.0  |
| H   | 1.67 | 1.1  | 1.24 | 1.94  | 1.0 | 0.45  |
| O   | 0.83 | 0.16  | 0.43 | 0.24  | 2.0 | 0.41  |

I have incorporated Tables 2, 3, 4 and 5 into the SOMnmR package as: "Terr_Nelson", "Aqua_Nelson", "Terr_Baldock" and "Aqua_Baldock", respectively.

[Back to the Main page](README.md)

