from classifiers import *
from subscript import *
from imblearn.under_sampling import RandomUnderSampler
from sklearn.model_selection import train_test_split

output_pth ='/User/predict_output/'


# use RRC(or RRPC) feature for prediction:
X_data, y_label = read_tsv('RRC_10.txt')
X, ind_X, y, ind_y = train_test_split(X_data, y_label, test_size=0.2, random_state=0)
rus = RandomUnderSampler(random_state=0)
X, y = rus.fit_resample(X, y)
print(Counter(y))
# imbalance method for testing set:
rus_indep = RandomUnderSampler(random_state=0, sampling_strategy=0.2)
ind_X, ind_y = rus_indep.fit_resample(ind_X, ind_y)
print(Counter(ind_y))

independent = np.zeros((ind_X.shape[0], ind_X.shape[1] + 1))
independent[:, 0], independent[:, 1:] = ind_y, ind_X
# SVM classifier
para_info, cv_res, ind_res = SVM_Classifier(X, y, indep=independent, fold=5)
# RF classifier
# para_info, cv_res, ind_res = RF_Classifier(X, y, indep=independent, fold=5)
# LR classifier
# para_info, cv_res, ind_res = LR_Classifier(X, y, indep=independent, fold=5)

# saving result for testing set
plot_roc_ind(independent, os.path.join(output_pth, 'ROC_ind.png'), label_column=0, score_column=2)
ind_metrics = calculate_metrics(independent[:, 0], independent[:, 2])
save_prediction_metrics_ind(ind_metrics, os.path.join(output_pth, 'metrics_IND.txt'))
