# coding=utf-8
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

'''Main function for UCI letter and spam datasets.
'''

# Necessary packages
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import argparse
import numpy as np
import struct
import time

from data_loader import data_loader
from gain import gain
from utils import rmse_loss


def main(args):
    '''Main function for UCI letter and spam datasets.

    Args:
      - data_name: letter or spam
      - miss_rate: probability of missing components
      - batch:size: batch size
      - hint_rate: hint rate
      - alpha: hyperparameter
      - iterations: iterations

    Returns:
      - imputed_data_x: imputed data
      - rmse: Root Mean Squared Error
    '''

    gain_parameters = {'batch_size': args.batch_size,
                       'hint_rate': args.hint_rate,
                       'alpha': args.alpha,
                       'iterations': args.iterations}

    #############################################
    # Load data and introduce missingness
    ##############################################
    daty_file = open("data/daty_Travel_Column.bin", 'rb')
    datr_file = open("data/datr_Travel_Column.bin", 'rb')
    ori_file = open("data/full_daty_column_binary.bin", 'rb')
    FHDI_file = open("data/final_daty_binary.bin", 'rb')
    nrow = 23772
    ncol = 50
    daty_list = []
    datr_list = []
    ori_list = []
    FHDI_list = []

    # -------------------------------
    # Read daty from a binary file
    # ------------------------------
    while True:
        data1 = daty_file.read(8)
        if not data1: break
        daty_list.append(struct.unpack('d', data1))

    daty_array = np.array(daty_list)
    daty = np.reshape(daty_array, (nrow, ncol), order='F')
    # np.savetxt('daty_python.txt', daty, fmt='%1.3f')

    # -----------------------------
    # Read datr from a binary file
    # -----------------------------
    while True:
        data2 = datr_file.read(4)
        if not data2: break
        datr_list.append(struct.unpack('i', data2))

    datr_array = np.array(datr_list)
    datr = np.reshape(datr_array, (nrow, ncol), order='F')
    # np.savetxt('datr_python.txt', datr, fmt='%d')

    # ------------------------------------------------
    # Read original complete data from a binary file
    # ------------------------------------------------
    while True:
        data3 = ori_file.read(8)
        if not data3: break
        ori_list.append(struct.unpack('d', data3))

    ori_array = np.array(ori_list)
    ori_data_x = np.reshape(ori_array, (nrow, ncol), order='F')
    np.savetxt('ori_data_x.txt', ori_data_x, fmt='%1.3f')

    # ------------------------------------------------
    # Read UP-FHDI file
    # ------------------------------------------------

    while True:
        data4 = FHDI_file.read(8)
        if not data4: break
        FHDI_list.append(struct.unpack('d', data4))

    FHDI_array = np.array(FHDI_list)
    FHDI_data_x = np.reshape(FHDI_array, (nrow, ncol))
    np.savetxt('FHDI_data_x.txt', FHDI_data_x, fmt='%1.3f')


    miss_data_x = daty.copy()
    miss_data_x[datr == 0] = np.nan
    np.savetxt('miss_data_x.txt', miss_data_x, fmt='%s')

    #####################################################

    # Impute missing data
    imputed_data_x = gain(miss_data_x, gain_parameters)

    daty_file.close()
    datr_file.close()
    ori_file.close()
    FHDI_file.close()


    # Report the RMSE performance
    rmse = rmse_loss(ori_data_x, FHDI_data_x, nrow, ncol)

    print()
    print('RMSE Performance: --------' + str(np.round(rmse, 6)) + '----------')
    print()


    return imputed_data_x


if __name__ == '__main__':
    start_time = time.time()
    # Inputs for the main function
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--data_name',
        choices=['letter', 'spam'],
        default='spam',
        type=str)
    parser.add_argument(
        '--miss_rate',
        help='missing data probability',
        default=0.2,
        type=float)
    parser.add_argument(
        '--batch_size',
        help='the number of samples in mini-batch',
        default=128,
        type=int)
    parser.add_argument(
        '--hint_rate',
        help='hint probability',
        default=0.9,
        type=float)
    parser.add_argument(
        '--alpha',
        help='hyperparameter',
        default=100,
        type=float)
    parser.add_argument(
        '--iterations',
        help='number of training interations',
        default=10000,
        type=int)

    args = parser.parse_args()

    # Calls main function
    imputed_data = main(args)
    # imputed_file = open("./imputed.txt", 'w')
    print('type of imputed is ', type(imputed_data))
    # print(imputed_data)
    np.savetxt('imputed.txt', imputed_data, fmt='%1.5f')
    print("------------------------ %s seconds --------------------" % (time.time() - start_time))
