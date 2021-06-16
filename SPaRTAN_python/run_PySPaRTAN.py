import sys
import os
import random
import numpy as np
import pandas as pd
from sklearn.model_selection import KFold
from sklearn.preprocessing import normalize
from scipy.io import loadmat
from scipy import spatial
from scipy import stats
from PySPaRTAN import PySPaRTAN

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--dataset_D", help="name of the dataset D which will be passed in", type=str, default="Dpbmc")
parser.add_argument("--dataset_P", help="name of the dataset P which will be passed in", type=str, default="Ppbmc5kn_CD8")
parser.add_argument("--dataset_Y", help="name of the dataset Y which will be passed in", type=str, default="Ypbmc5kn_CD8")
parser.add_argument("--spectrumP", help="Dimension reduction coefficient on protein space", type=float, default=0.7)
parser.add_argument("--lamda", help="LASSO regression coefficient", type=float, default=0.001)
parser.add_argument("--rsL2", help="ridge regression coefficient", type=float, default=0.001)
parser.add_argument("--normalization", help="type of normalizion, no normalization if set to empty", type=str, default="l2")
parser.add_argument("--input_dir", help="directory of input files", type=str, default="../data/inputs")
parser.add_argument("--output_dir", help="directory of output files", type=str, default="../data/outputs")
parser.add_argument('--fold', help="how many folds for the cross_validation.No cross_validation and using default/specified parameters if set to 0", type=int, default=0)

args = parser.parse_args()



if not os.path.exists(args.output_dir):
    os.makedirs(args.output_dir)

print("Read in datasets D, P, and Y ...")
D_ori = pd.read_csv(os.path.join(args.input_dir, args.dataset_D+'.csv'), index_col=0)
P_ori = pd.read_csv(os.path.join(args.input_dir, args.dataset_P+'.csv'), index_col=0)
Y_ori = pd.read_csv(os.path.join(args.input_dir, args.dataset_Y+'.csv'), index_col=0)

TF_name = list(D_ori.columns)
cell_name = list(Y_ori.columns)
gene_name = list(Y_ori.index)
protein_name = list(P_ori.columns)

D_mat = D_ori.values
P_mat = P_ori.values
Y_mat = Y_ori.values

# normalize the dataset
if args.normalization != "":
    D = normalize(D_mat, norm=args.normalization, axis=0)
    Y = normalize(Y_mat, norm=args.normalization, axis=0)
    P = normalize(P_mat, norm=args.normalization, axis=1)

# create the object of SPaRTAN
reg = PySPaRTAN()



# cross-validate to determine hyperparameters
fold = args.fold
if fold != 0: #using cross validation to determine the optimal parameters
    
    lamdas = [0.001, 0.01, 0.1, 0.2, 0.3]
    rsL2s = [0, 0.001, 0.01]
    spectrumAs = [1]
    spectrumBs = [0.5, 0.6, 0.7]

    lenlamdas = len(lamdas)
    lenrsL2s = len(rsL2s)
    lenspAs = len(spectrumAs)
    lenspBs = len(spectrumBs)
    
    corr_all_spearman = np.zeros((lenspAs, lenspBs, lenlamdas, lenrsL2s))
    for a in range(0, lenspAs):
        for b in range(0, lenspBs):
            for l in range(0, lenlamdas):
                for r in range(0, lenrsL2s):
                    print(
                        "cross validating spectrumA={}, spectrumB={}, lambda={}, rsL2={}".format(
                            spectrumAs[a], spectrumBs[b], lamdas[l], rsL2s[r]
                        )
                    )
                    sum_corr_spearman = 0
    
                    kf = KFold(n_splits=fold)
                    for train_index, test_index in kf.split(P_mat):
    
                        # split the data into train and test set
                        P_train, P_test = P_mat[train_index, :], P_mat[test_index, :]
                        Y_train, Y_test = Y_mat[:, train_index], Y_mat[:, test_index]
    
                        # normalize the train and test set
                        Y_train = normalize(Y_train, axis=0)
                        Y_test = normalize(Y_test, axis=0)
    
                        P_train = normalize(P_train, axis=1)
                        P_test = normalize(P_test, axis=1)
    
                        # train the model
                        reg.fit(
                            D,
                            P_train,
                            Y_train,
                            lamda=lamdas[l],
                            rsL2=rsL2s[r],
                            spectrumA=spectrumAs[a],
                            spectrumB=spectrumBs[b],
                        )
    
                        # get predicted value Y_pred  on P_test
                        Y_pred = reg.predict(P_test)
    
                        # get the correlation bewteen Y_pred and Y_test
                        corr_spearman = reg.get_corr(Y_pred, Y_test)
    
                        sum_corr_spearman = sum_corr_spearman + corr_spearman
    
                    corr_all_spearman[a, b, l, r] = sum_corr_spearman / fold
                    
   
    # retrive the best parameters
    max_a, max_b, max_l, max_r = np.unravel_index(
        corr_all_spearman.argmax(), corr_all_spearman.shape
    )
    
    lamda_best = lamdas[max_l]
    rsL2_best = rsL2s[max_r]
    spectrumA_best = spectrumAs[max_a]
    spectrumB_best = spectrumBs[max_b]
       
    print("lamda_best={}, rsL2_best={}, spectrumA_best={},spectrumB_best={}".format(lamda_best,rsL2_best,spectrumA_best,spectrumB_best))
             
                    
else: #fold ==0: using default/specified paramters
 
    lamda_best=args.lamda
    rsL2_best=args.rsL2
    spectrumA_best=1
    spectrumB_best=args.spectrumP

print("Processing ...")
# re-train the model
reg.fit(D, P, Y, lamda_best, rsL2_best, spectrumA_best, spectrumB_best)

# retrieve W, projD, projP
W = reg.get_W()
projD = reg.get_projD()
projP = reg.get_projP()


df_W = pd.DataFrame(data=W, index=TF_name, columns=protein_name)
df_projP = pd.DataFrame(data=projP, index=cell_name, columns=protein_name)
df_projD = pd.DataFrame(data=projD, index=TF_name, columns=cell_name)

outfile_W = os.path.join(args.output_dir, "W.csv")
outfile_projP = os.path.join(args.output_dir, "projP.csv")
outfile_projD = os.path.join(args.output_dir, "projD.csv")

df_W.to_csv(outfile_W)
df_projP.to_csv(outfile_projP)
df_projD.to_csv(outfile_projD)

# with open(outfile_W, "w") as outfile:
    # np.savetxt(outfile, W, fmt="%-7.4f")

# with open(outfile_projP, "w") as outfile:
    # np.savetxt(outfile, projP, fmt="%-7.4f")

# with open(outfile_projD, "w") as outfile:
    # np.savetxt(outfile, projD, fmt="%-7.4f")

print("Process finished successfully!")
