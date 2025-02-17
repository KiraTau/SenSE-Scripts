# SenSE-Scripts
This repository holds scripts developed by Chukwuemeka Ike as part of his circadian estimation research. The folders are divided into two distinct types:
* Algorithms - scripts in these folders run the major algorithms
* Support - combination of utility scripts, Simulink models, and data files that are commonly used across the algorithms

Note: The data need to run the example scripts are all located in a Data/ folder which Zidi Tao has access to. All requests for access can be directed to either him or Dr. Agung Julius.


## Algorithms
### AutoEncoder + Gated Recurrent Unit (*AE_GRU/*)
Implementation of the black-box identification work done with Yunshi Wen in [Personalized Data-Driven State Models of the Circadian Dynamics in a Biometric Signal](). The system uses an autoencoder to compress a history of actigraphy and light exposure into a lower-dimensional space, then uses a gated recurrent unit as a state transition model within that latent space.

#### Citation
```bibtex
@inproceedings{ike2024personalized,
  title={Personalized Data-Driven State Models of the Circadian Dynamics in a Biometric Signal},
  author={Ike, Chukwuemeka O and Wen, Yunshi and Wen, John T and Oishi, Meeko M K and Brown, Lee K and Julius, A Agung},
  booktitle={2024 46th Annual International Conference of the IEEE Engineering in Medicine \& Biology Society (EMBC)},
  pages={1--5},
  year={2024},
    month={7},
}
```

### Kalman Filter (*KF/*)
These scripts implement the circadian phase shift estimation algorithm presented in [Efficient Estimation of the Human Circadian Phase via Kalman Filtering](https://scholar.google.com/citations?view_op=view_citation&hl=en&user=kgLMmEIAAAAJ&citation_for_view=kgLMmEIAAAAJ:UeHWp8X0CEIC). 

#### Citation
```bibtex
@inproceedings{ike2023efficient,
  title={Efficient Estimation of the Human Circadian Phase via Kalman Filtering},
  author={Ike, Chukwuemeka O and Wen, John T and Oishi, Meeko M K and Brown, Lee K and Julius, A Agung},
  booktitle={2023 45th Annual International Conference of the IEEE Engineering in Medicine \& Biology Society (EMBC)},
  pages={1--6},
  year={2023},
    month={7},
}
```
### Observer-based Filter (*OBF/*)
Theh observer-based filter was developed for circadian phase shift estimation using accessible biometric signals in [Fast tuning of observer-based circadian phase estimator using biometric data](https://scholar.google.com/citations?view_op=view_citation&hl=en&user=kgLMmEIAAAAJ&citation_for_view=kgLMmEIAAAAJ:d1gkVwhDpl0C).

#### Citation
``` bibtex
@article{ike2022fast,
  title={Fast tuning of observer-based circadian phase estimator using biometric data},
  author={Ike, Chukwuemeka O and Wen, John T and Oishi, Meeko MK and Brown, Lee K and Julius, A Agung},
  journal={Heliyon},
  volume={8},
  number={12},
  year={2022},
  publisher={Elsevier}
}
```

## Particle Filter (*PF/*)
The particle filter was used in building a circadian state estimation framework.

#### Citation
``` bibtex
@article{ike2024model,
  title={Model-Based Human Circadian State Estimation with Wearable Device Data},
  author={Ike, Chukwuemeka O and Wen, John T and Oishi, Meeko MK and Brown, Lee K and Julius, A Agung},
}
```

## Synthetic Data Generation (*Synth/*)
Scripts for generating synthetic circadian heart rate data. The process of generating a set of signals is included in the folder's README. To cite, use the citation for the Particle Filter paper.




# Support Code
## Models
This folder contains Simulink and MATLAB implementations of models of:
* the electrical activity in the heart,
* the circadian system, and
* the OBF and the ANF (Yin 2020)

Each script cites the paper it is based on and any modifications made to the original work. They're used all over the algorithm scripts for several purposes.


## Utils
This folder contains utility scripts for such things as:
* loading biometric data
* calculating filter-based phase shifts
* converting activity in steps to light based on the method by Huang et al. (2019)
* generating a simulated light signal




# Fitbit Python
The lone script in this folder allows a user to download their own Fitbit data over a specified length of time. The script contains commented instructions on how to edit it to get what you need.