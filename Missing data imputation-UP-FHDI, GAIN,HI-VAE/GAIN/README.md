# Codebase for "Generative Adversarial Imputation Networks (GAIN)"

Authors: Jinsung Yoon, James Jordon, Mihaela van der Schaar

Paper: Jinsung Yoon, James Jordon, Mihaela van der Schaar, 
"GAIN: Missing Data Imputation using Generative Adversarial Nets," 
International Conference on Machine Learning (ICML), 2018.
 
Paper Link: http://proceedings.mlr.press/v80/yoon18a/yoon18a.pdf

Contact: jsyoon0823@gmail.com

This directory contains implementations of GAIN framework for imputation
using two UCI datasets.

-   UCI Letter (https://archive.ics.uci.edu/ml/datasets/Letter+Recognition)
-   UCI Spam (https://archive.ics.uci.edu/ml/datasets/Spambase)

To run the pipeline for training and evaluation on GAIN framwork, simply run 
python3 -m main_letter_spam.py.

Note that any model architecture can be used as the generator and 
discriminator model such as multi-layer perceptrons or CNNs. 

### Command inputs:

-   data_name: letter or spam
-   miss_rate: probability of missing components
-   batch_size: batch size
-   hint_rate: hint rate
-   alpha: hyperparameter
-   iterations: iterations

### Example command

```shell
$ python3 main_letter_spam.py --data_name spam  miss_rate 0.2 --batch_size 128 --hint_rate 0.9 --alpha 100 iterations 10000
```

### Outputs

-   imputed_data_x: imputed data
-   rmse: Root Mean Squared Error

### GAIN for incomplete binary data
Original GAIN is adapted from Codebase for "Generative Adversarial Imputation Networks (GAIN)" on Github released by the authors (Jinsung Yoon, James Jordon, Mihaela van der Schaar)
To run modified GAIN for an incomplete binary data, simply run

$ python3 main_letter_spam.py --iterations 1000

Note that default parameters are used except for the number of iterations. iterations can be reduced for a tolerable execution time. The input data dimensions and directory must
be specified in main_letter_spam.py