"""
This script intends to use PySPaRTAN module to generate predicted matrices used in the paper.

PySPaRTAN has 3 Hyperparameters:Spectrum P, lamda, and rsL2
We can run PySPaRTAN by specifying some values to those parameters or using default ones in the script.
We can also use cross-validation to generate the optional values for those
Hyperparameters at first, and then run PySPaRTAN to generate the projections.


When running this script from command line, the following parameters can be added to the command:
    --input_dir : directory of input files, default="../data/inputs"
    --output_dir : directory of output files, default="../data/outputs"
    --dataset_D : the name of gene-TF matrix file.
                  the file requires .csv format.
                  only contains name here, do not include ".csv" extension
    --dataset_P : the name of protein matrix
                  the same other requirements as --dataset_D
    --dataset_Y : the name of gene expression matrix
                  the same other requirements as --dataset_D
    --spectrumP : Dimension reduction coefficient on protein space, default=0.7  
    --lamda :  LASSO regression coefficient,default=0.001
    --rsL2 : ridge regression coefficient, default=0.001
    --normalization : type of normalizion performed on matrices,
                      default is l2 normalization
    --fold : how many folds to be used when doing cross-validation.
             default=0, means using specified hyper-parameters, do not do cross-validation
             
             
This script requires numpy, pandas, sklearn to be installed in the python running environment

"""
import argparse
import os
import numpy as np
import pandas as pd
from sklearn.model_selection import KFold
from sklearn.preprocessing import normalize
from PySPaRTAN import PySPaRTAN

parser = argparse.ArgumentParser()

parser.add_argument("--input_dir", help="directory of input files",
                    type=str, default="../data/inputs")
parser.add_argument("--output_dir", help="directory of output files",
                    type=str, default="../data/outputs")
parser.add_argument("--dataset_D", help="name of gene-TF matrix",
                    type=str, default="Dpbmc")
parser.add_argument("--dataset_P", help="name of the dataset P which will be passed in",
                    type=str, default="Ppbmc5kn_CD16")
parser.add_argument("--dataset_Y", help="name of the dataset Y which will be passed in",
                    type=str, default="Ypbmc5kn_CD16")
parser.add_argument("--spectrumP", help="Dimension reduction coefficient on protein space",
                    type=float, default=0.7)
parser.add_argument("--lamda", help="LASSO regression coefficient",
                    type=float, default=0.001)
parser.add_argument("--rsL2", help="ridge regression coefficient",
                    type=float, default=0.001)
parser.add_argument("--normalization", help="type of normalizion performed on matrices,\
                    no normalization if set to empty", type=str, default="l2")
parser.add_argument('--fold', help="how many folds for the cross_validation.\
                    No cross_validation and using default/specified parameters if set to 0",
                    type=int, default=5)

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

# cross-validate to determine optimal parameters
fold = args.fold
if fold != 0:  # using cross validation to determine the optimal parameters
    
    lamdas = [0.001, 0.01, 0.1, 0.2, 0.3]
    rsL2s = [0.001, 0.01, 0.1]
    spectrumAs = [1]
    spectrumBs = [0.5, 0.6, 0.7]

    lenlamdas = len(lamdas)
    lenrsL2s = len(rsL2s)
    lenspAs = len(spectrumAs)
    lenspBs = len(spectrumBs)
    
    corr_all_pearson = np.zeros((lenspAs, lenspBs, lenlamdas, lenrsL2s))
    for a in range(0, lenspAs):
        for b in range(0, lenspBs):
            for l in range(0, lenlamdas):
                for r in range(0, lenrsL2s):
                    print("cross validating spectrumP={}, lambda={}, rsL2={}"
                          .format(spectrumBs[b], lamdas[l], rsL2s[r]))
                    sum_corr_pearson = 0
    
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
                        corr_pearson = reg.get_corr(Y_pred, Y_test)
    
                        sum_corr_pearson = sum_corr_pearson + corr_pearson
    
                    corr_all_pearson[a, b, l, r] = sum_corr_pearson / fold

    # retrive the best parameters
    max_a, max_b, max_l, max_r = np.unravel_index(
        corr_all_pearson.argmax(), corr_all_pearson.shape
    )
    
    lamda_best = lamdas[max_l]
    rsL2_best = rsL2s[max_r]
    spectrumA_best = spectrumAs[max_a]
    spectrumB_best = spectrumBs[max_b]
       
    print("lamda_best={}, rsL2_best={}, spectrumA_best={}, spectrumB_best={}"
          .format(lamda_best, rsL2_best, spectrumA_best, spectrumB_best))

else:  # fold ==0: using default/specified paramters

    lamda_best = args.lamda
    rsL2_best = args.rsL2
    spectrumA_best = 1
    spectrumB_best = args.spectrumP

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


print("Process finished successfully!")
