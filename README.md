# **Data preprocessing-Create missingness**

Let $\textbf{U}(n,p)$ denote a complete real-world dataset, and its variables are indexed by $0,1,\dots,p-1$. In line 2 of Algorithm, split $\textbf{U}(n, p)$ as disjoint $\textbf{U}_1(0.7n, p)$ and $\textbf{U}_2(0.3n, p)$ such that $\textbf{U}_1 \cup \textbf{U}_2 = \textbf{U}$. 

Create missingness to $\textbf{U}_1(0.7n, p)$ with different missing mechanisms, i.e., missing completely at random (MCAR) or missing not at random (MNAR)

## MCAR

The MCAR mechanism indicates that the occurrence of missing data is completely independent of the observed data; in other words, the missing data are completely random.

### Usage

#### Command

For linux

```linux
1. source /opt/intel/oneapi/setvars.sh
```

```linux
2. mpiicc -std=c++11 -o main_MPI Unbiased_Missingness.cpp
```

```linux
3. sbatch run.sbatch
```

## MNAR

The MNAR mechanism is more complex, indicating that the occurrence of missing data is related not only to the observed data but also to the values of the missing data themselves. In real-world scenarios, data missingness often follows the MNAR mechanism.

### Usage

#### Command

For linux

```linux
1. source /opt/intel/oneapi/setvars.sh
```

```linux
2. mpiicc -std=c++11 -o main_MPI Biased_Missingness.cpp
```

```linux
3. sbatch run.sbatch
```

# **Missing data imputation-UP-FHDI, GAIN,HI-VAE**

Cure incomplete data $\textbf{U}_1(0.7n, p)$ using different imputation techniques to obtain the cured data $\widehat{\textbf{U}}_1(0.7n, p)$.

## **UP-FHDI**

Ultra data-oriented parallel fractional hot-deck imputation (UP-FHDI) is a general-purpose, assumption-free imputation software capable of curing big incomplete data with complex and irregular missing patterns. It inherits all the strengths of the _R_ Package FHDI. UP-FHDI can tackle datasets with up to one million instances, 10,000 variables, and 30% missing values.

- [Benchmarks](#Benchmarks)
- [Usage](#Usage)
  - [Dependencies](#Dependencies)
  - [Command](#Command)
- [Citation](#Citation)
- [Acknowledgements](#Acknowledgements)



Please see a tutorial video in [UP-FHDI on HPC](https://www.youtube.com/watch?v=Dr0x5lZsVuU) to illustrate the use of UP-FHDI with example datasets.

### Usage

#### Dependencies

- Intel MPI
- Access to HPC Facilities


#### Command
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

## **GAIN**

GAIN (Generative Adversarial Imputation Nets) is a deep learning-based method designed to impute missing data by leveraging the Generative Adversarial Networks (GAN) framework. Introduced by Jinsung Yoon, James Jordon, and Mihaela van der Schaar in 2018, GAIN adapts GANs to the context of data imputation/ref.

### Usage

The detailed usage can be found in the README.md file in this folder.


## **HI-VAE**

Handling Incomplete Heterogeneous Data using VAEs refers to the use of Variational Autoencoders (VAEs) to model and impute datasets that contain both missing values and different types of variables, such as continuous, categorical, ordinal, and count data. This approach, often exemplified by the HI-VAE model, extends the standard VAE framework to support mixed data types and naturally handle missing entries during training. By learning a shared latent representation, VAEs can generate plausible imputations and improve downstream tasks like classification or clustering on incomplete and diverse datasets.

### Usage

The detailed usage can be found in the README.md file in this folder.


# **Two-staged feature selection**

Growing data volume will result in substantial training time for ML models. Dimension reduction is an essential technique to ease the model training process. Suppose we need to have a subset of $w$ features from $\textbf{U}_1(0.7n, p)$, $\widehat{\textbf{U}}_1(0.7n, p)$, and $\textbf{U}_2(0.3n, p)$ such that $w \ll p$. The IPIF utilizes the two-step feature selection method proposed in \cite{Yang:2023} for effective variable reduction in line 10 of Algorithm. 

The proposed feature selection method incorporates predictor-target correlations and inter-predictor correlations in a two-step process.Let $\textbf{y}_I = \{y_0,\dots, y_{p-2}\}$ as predictors and $y_{p-1}$ as the target variable. Discretize $\textbf{y}_I$ to $\textbf{z}_I$ by categories $\boldsymbol{G}$ using the estimated sample quantiles. Similarly, discrete $y_{p-1}$ to $z_{p-1}$ by the category $G_{p-1}$. The first step filters $v_I$ important features based on mutual information by:

\begin{flalign}
      \boldsymbol{MI} = \sum_{g=1}^{G_l} \sum_{k=1}^{G_{p-1}} p(g,k)  \log \frac{p(g,k)}{p(g)p(k)} 
  \end{flalign}
  
where $G_l \in \boldsymbol{G}$, $p(g)$ and $p(k)$ are marginal probability, and $p(g,k)$ is joint probability. $v_I$ features with top mutual information with respect to $\boldsymbol{MI}$ will be selected such that $v_I < p$. 

The second step leverages graphical models based on the inverse covariance matrix (denoted as $\boldsymbol{\Sigma}^{-1}$) to further reduce $v_I$ features to a final small-sized subset of $w$ features. $\boldsymbol{\Sigma}^{-1}$ displays the partial correlations of variables, and zero elements imply the conditional independency given the rest. We adopt the R package $glasso$ for the estimation of a sparse $\boldsymbol{\Sigma}^{-1}$ using a lasso penalty. Note that users may skip this stage when data volume is moderate.

## Parallel MI


## Glasso

# **Downstream predicitve models**

The downstream predictive models include MLP, SVM, and Random Forest, which were applied to the SWARM, FMA, and COVID datasets, respectively. Except for the MLP model, which was configured with four hidden layers of sizes (50, 50, 50, 25), all other model parameters followed the default settings provided by the Python scikit-learn library.


# **Dataset**
| Dataset  | # Instances | # Variables |   Source |
| :---: | :---: | :---: | :---: |
| Swarm | 24016 | 2401 | UCI |
| FMA | 106574 | 518 | UCI |
| Covid | 1056661 | 768 | Kaggle |

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

@Article{Yoon:2018,
  author = {Jinsung Yoon and James Jordon and Mihaela van der Schaar},
  title = {GAIN: Missing data imputation using generative adversarial nets},
  journal = {35th International Conference on Machine Learning},
  year = {2018},
  pages = {5689--5698},
}

@Article{HI_VAE:2020,
  author = {Alfredo Nazabal and Pablo M. Olmos and Zoubin Ghahramani and Isabel Valera},
  title = {Handling Incomplete Heterogeneous Data using VAEs},
  journal = {Pattern Recognition},
  volume={107},
  year={2020},
  number={107501},
}

```
