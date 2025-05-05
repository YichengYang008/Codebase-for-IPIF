# **UP-FHDI**

Ultra data-oriented parallel fractional hot-deck imputation (UP-FHDI) is a general-purpose, assumption-free imputation software capable of curing big incomplete data with complex and irregular missing patterns. It inherits all the strengths of the _R_ Package FHDI. UP-FHDI can tackle datasets with up to one million instances, 10,000 variables, and 30% missing values.

- [Benchmarks](#Benchmarks)
- [Usage](#Usage)
  - [Dependencies](#Dependencies)
  - [Command](#Command)
- [Citation](#Citation)
- [Acknowledgements](#Acknowledgements)



Please see a tutorial video in [UP-FHDI on HPC](https://www.youtube.com/watch?v=Dr0x5lZsVuU) to illustrate the use of UP-FHDI with example datasets.

# Benchmarks
Benchmarks with missing values are publicly accessible in [IEEE DataPort](https://ieee-dataport.org/open-access/incomplete-datasets-fhdi).

| Dataset  | # Instances | # Variables | Missing rate |   Source |
| :---: | :---: | :---: | :---: |  :---: |
| Synthetic 1  | 100  | 80 | 30% | Synthetic |
| Synthetic 2  | 10000  | 0.1M | 30% |  Synthetic |
| Synthetic 3 | 0.1M  | 10000 | 30% | Synthetic |
| Earthquake | 901512  | 15 | 30% | USGS |
| Bridge Strain | 492641  | 31 | 30% | InTrans |
| Travel Time | 23772  | 50 | 30% | IEEE DataPort |
| CT Slices | 53500  | 380 | 30% | UCI |
| Swarm | 24016 | 2400 | 30% | UCI |
| p53 | 31159 | 5408 | 30% | UCI |
| Radar | 325834 | 175 | 30% | UCI |

# Usage

## Dependencies

- Intel MPI
- Access to HPC Facilities


## Command
UP-FHDI has been broadly validated on two high-performance computing (HPC) facilities: Condo2017 and Stampede2. Condo2017 is a free HPC platform open to Iowa State University faculties. Users can run the following commands to deploy UP-FHDI on Condo2017. 

Load the Intel compiler:
```linux
module load intel/18.3
```

Compile UP-FHDI program:
```linux
mpiicc –o main_MPI main_MPI.cpp
```

Submit a job script to launch MPI code with requested computing resources
```linux
sbatch run.sbatch
```

Stampede2 is a flagship HPC facility at the University of Texas at Austin open to global researchers with awarded allocations. Users can run the following commands to deploy UP-FHDI on Stampede2. 

Compile UP-FHDI program:
```linux
mpicxx –o main_MPI main_MPI.cpp
```

Submit a job script to launch MPI code with requested computing resources
```linux
sbatch run.sbatch
```

Note that users should carefully investigate how to deploy UP-FHDI on other HPC facilities! 

# Citation
Please kindly cite the following papers if you find this software useful. Thanks!
- Yicheng Yang, Yonghyun Kwon, Jaekwang Kim, and In Ho Cho, 2023. [Ultra data-oriented parallel fractional hot-deck imputation with efficient linearized variance estimation](https://ieeexplore.ieee.org/document/10054160),  _IEEE Transactions on Knowledge and Data Engineering_ (accepted).
- Yicheng Yang, Jaekwang Kim, and In Ho Cho, 2020. [Parallel fractional hot-deck imputation and variance estimation for big incomplete data curing](https://ieeexplore.ieee.org/document/9214981), _IEEE Transactions on Knowledge and Data Engineering_ 34(8), 3912-3926 [DOI: 10.1109/TKDE.2020.3029146].
- Jongho Im,	In Ho Cho, and Jaekwang Kim, 2018. [FHDI: an R package for fractional hot-deck imputation for multivariate missing data](https://journal.r-project.org/archive/2018/RJ-2018-020/index.html), _The R Journal_ 10(1), 140-154 [DOI: 10.32614/RJ-2018-020].

```latex
@Article{Yicheng:2023, 
  author = {Yicheng Yang and Yonghyun Kwon and Jaekwang Kim and In Ho Cho},
  title = {Ultra data-oriented parallel fractional hot-deck
  imputation with efficient linearized variance estimation},
  journal = {IEEE Transactions on Knowledge and Data Engineering (accepted)},
  year = {2023},
}

@Article{Yicheng:2020,
  author = {Yicheng Yang and Jaekwang Kim and In Ho Cho},
  title = {Parallel fractional hot deck imputation and variance estimation for big incomplete data curing},
  journal = {IEEE Transactions on Knowledge and Data Engineering},
  year = {2020},
  volume = {34},
  number = {8},
  pages = {3912--3926},
  doi = {10.1109/TKDE.2020.3029146},
}

@Article{RJ-2018-020,
  author = {Jongho Im and In Ho Cho and Jae Kwang Kim},
  title = {FHDI: An R Package for Fractional Hot Deck Imputation},
  year = {2018},
  journal = {The R Journal},
  url = {https://doi.org/10.32614/RJ-2018-020},
  volume = {10},
  number = {1},
  pages = {140--154},
  doi = {10.32614/RJ-2018-020},
}
```

# Acknowledgements

This software is supported by National Science Foundation
(NSF) grant number OAC-1931380. The high-performance
computing facility used for this research is partially supported
by the HPC@ISU equipment at ISU, some of which have been
purchased through funding provided by NSF CNS 1229081
and CRI 1205413. Ultra data applications of this software used
the Extreme Science and Engineering Discovery Environment
(XSEDE), NSF ACI-1548562.