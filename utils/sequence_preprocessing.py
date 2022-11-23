#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
@File: sequence_preprocessing.py
@Time: 2022/5/2 22:41
@Author: genqiang_wu@163.com
@desc: 

"""

def get_sequence_samples(Seq, window_size):
    """
    处理单条序列数据,将单条数据进行切片，返回切片序列和修饰位点在序列中的索引
    """
    aminodata = str(Seq)
    left_window_size = window_size // 2
    sequences = []
    k_index = []
    for r in range(len(aminodata)):
        if aminodata[r] == 'K':
            k_index.append(r)
            if len(aminodata[r+1:]) < left_window_size:  # 判断右边长度是否小于left_window_size
                aminovecR = aminodata[r:] + 'X' * (left_window_size - len(aminodata[r + 1:]))
                aminovecL = aminodata[r - left_window_size:r]
            elif len(aminodata[:r]) < left_window_size:  # 判断前面长度是否小于left_window_size
                aminovecL = 'X' * (left_window_size - len(aminodata[:r])) + aminodata[:r]
                aminovecR = aminodata[r:r + left_window_size + 1]
            else:  # 正常情况
                aminovecL = aminodata[r - left_window_size:r]
                aminovecR = aminodata[r:r + left_window_size + 1]
            aminovec = aminovecL + aminovecR
            sequences.append(aminovec)
    return sequences, k_index # sequences: 切分好的序列. k_index: K 的索引值

