#!/usr/bin/python
# -*- coding: utf-8 -*-
# 対称行列のMatrixMerket形式ファイルから入力ファイルを作成

import sys
args=sys.argv
argc=len(args) # 引数の数

if  argc == 3 :
#    print(args)  
    print("input  file :: " + args[1]) # 第一引数 Matrix Market 入力ファイル
    print("output file :: " + args[2]) # 第二引数 COO Matrix    出力ファイル
else:
    exit()
    
# print("START Matrix Market --> CRS Matrix ") 
mmfile  = open(args[1],"r")
coofile = open(args[2],"w+")

lines=mmfile.readlines()
mmlines=[]
for line in lines:
    if line[0] != "%":
        mmlines.append( list(line.split(" ")) )

#####
# 行,列,非零要素出力
# outline = mmlines[0][0]+" "+mmlines[0][1]+" "+mmlines[0][2]
# 対称性なので非ゼロ要素は 2*nnz-n  !! 列成分が全て零の行がある場合は不正
outline = mmlines[0][0]+" "+mmlines[0][1]+" "+str(2*int(mmlines[0][2]) - int( mmlines[0][0]) )+"\n"
coofile.write(outline)

#####
# coo形式の行列を出力
cooline=[]
for line in mmlines[1:len(mmlines)]:
    cooline.append([int(line[0]),int(line[1]),float(line[2])])
    if int(line[0]) != int(line[1]) :
        cooline.append([int(line[1]),int(line[0]),float(line[2])])

for line in sorted(cooline[0:len(cooline)]):
    outline = str(line[0])+" "+str(line[1])+" "+str(line[2])+"\n"
#    print outline
    coofile.write(outline)

#####

mmfile.close()
coofile.close()
