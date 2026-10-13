import numpy as np
from sklearn.model_selection import GridSearchCV
from group_lasso import GroupLasso
import zarr
import sys
from sklearn.linear_model import LassoCV
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.linear_model import Lasso

#load X and Y
X=sys.argv[1]
Y=sys.argv[2]
out_dir=sys.argv[3]
name=sys.argv[4]

X=zarr.load(X)
Y=zarr.load(Y)

alphas = np.logspace(-4, 1, 100)  # wider range

model = make_pipeline(
    StandardScaler(),
    LassoCV(alphas=alphas, cv=5, max_iter=1000)
)


print("Fitting model")
model.fit(X,Y)

lasso = model.named_steps['lassocv']
print("Best alpha:", lasso.alpha_)

alpha=lasso.alpha_

betas_lasso=np.zeros((X.shape[1],100))
for iteration in range(0,100):
    print ("Iteration: ",iteration)
    lasso = Lasso(alpha=alpha)
    lasso.fit(X,Y)

    beta=lasso.coef_
    betas_lasso[:,iteration]=beta
    print("Number of non 0 elements: " ,len(np.where(beta!=0)[0]))

mean_beta=np.mean(betas_lasso,axis=1).tolist()
pip=np.count_nonzero(betas_lasso, axis=1)/100

#save to file
zarr.save(out_dir+"betas_"+name+".zarr",betas_lasso)

with open(out_dir+"meanBeta_"+name+".txt", "a+") as f:
     s = ','.join(str(x) for x in mean_beta)
     f.write(s+"\n")

with open(out_dir+"pip_"+name+".txt", "a+") as f:
     s = ','.join(str(x) for x in pip)
     f.write(s+"\n")
