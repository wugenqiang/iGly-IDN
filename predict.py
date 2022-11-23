#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
@File: predict.py
@Time: 2022/5/3 20:41
@Author: genqiang_wu@163.com
@desc:

The model to predict the glycation sites.

"""

from keras.models import load_model
from Bio import SeqIO # pip install biopython # 导入SeqIO模块

from utils.feature_extraction import one_hot
from utils.sequence_preprocessing import get_sequence_samples

'''
处理example中fasta格式的数据并进行位点预测
步骤如下：
1.根据窗口大小切分待预测数据
2.特征提取
3.加载模型
4.得到预测值
5.用字典存储预测数据
6.将预测数据写入csv文件

'''

import numpy as np
import pandas as pd

window_size = 31 # 序列窗口大小设置为31
sequences_need_predict_fasta = 'example/sequences_need_to_be_predicted_example.fasta' # 待预测序列

pre_site = [] # 预测为糖基化位点的k的位置列表
count = 0 # 预测蛋白质序列条数

for seq_i in SeqIO.parse(sequences_need_predict_fasta, 'fasta'):
    # print(seq_i)
    # Step 1: 根据窗口大小切分数据, window size = 31
    sequences, k_index = get_sequence_samples(seq_i.seq, window_size)  # 对待预测序列进行切分
    # print(sequences, k_index)

    # Step 2: 特征提取
    sequence_feature = one_hot(sequences, window_size)  # one hot

    # Step 3: 加载模型
    model = load_model('models/iGly-IDN_model.h5')

    # Step 4: 得到预测值
    sequences_pred = model.predict(sequence_feature, verbose=1)
    # print(sequences_pred[:, 1])
    y_pred = sequences_pred[:, 1] # 预测概率值

    k_index_pred_true = []  # 预测为真实值的索引

    for i in range(len(y_pred)):
        if y_pred[i] >= 0.5:
            k_index_pred_true.append(k_index[i])

    num_k_index_pred_true = len(k_index_pred_true)  # 预测为真实值的个数
    # print(k_index_pred_true, num_k_index_pred_true)

    k_loc_pred_true = [loc+1 for loc in k_index_pred_true] # 预测为真实值的k的位置
    # print(k_loc_pred_true)

    # Step 5: 用字典存储预测数据
    content = {
        'number': count,
        'name': seq_i.name,
        'sequence': seq_i.seq,
        'pre_site': k_loc_pred_true,
    }
    count += 1
    pre_site.append(content)
    print(content)

# Step 6: 将预测数据写入csv文件
pd.DataFrame(pre_site).to_csv('result/predict_result/example_pre_site.csv', index=False)



