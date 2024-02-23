# StabilityCCA

This repository contains MATLAB implementations of two sparse CCA methods, **SCCA-EC** and **PMD-CCA**, and a method for performing stability selection using sparse CCA, **StabilityCCA**.

PMD-CCA is the penalised matrix decomposition approach presented in [1].
SCCA-EC combines the L_1 and L_2 -norm constraints from [1] with the alternating projected gradient algorithm from [2].
Both methods can be used to calculate regularisation paths: if no sparsity
parameters are given, the functions return a sequence of CCA models, each for a pair of sparsity parameter values. StabilityCCA uses this feature to
calculate stability paths (see [3] and [4]): data is subsampled in order to estimate the probability that each variable would be selected by SCCA at a given
level of sparsity.

The main functions are `SCCA-EC/SCCAec`, `PMD-CCA/PMDCCA` and `stabilityCCA`.

"Installing" a MATLAB package is easy: just clone this repository and add
it and its subfolders to your MATLAB path:

```
addpath(genpath('path_to/stabilityCCA'))
```

## SCCA-EC

Create some simulated data, run `SCCAec` and plot the results
```
create_stabilityCCA_example
load stabCCA_example

[A,B,~,~,~,cxy] = SCCAec(data2.X,data2.Y);
figure
    subplot(2,1,1);plotCCApath(A,cxy(:,1),find(data2.A));title('X-view')
    subplot(2,1,2);plotCCApath(B,cxy(:,2),find(data2.B));title('Y-view')
```
The resulting plots show regularisation paths for all variables from the two views with the ground truth variables highlighted.

If you want "normal" sparse CCA, provide values for sparsity parameters:
```
cx = 3; % X-view sparsity
cy = 4; %Y-view sparsity
[A,B,r,U,V] = SCCAec(data1.X,data1.Y,'cxy',[cx,cy]);

figure
    subplot(2,2,1);bar(A);
    xlabel('variable');ylabel('coefficient');title('X-view')
    subplot(2,2,2);bar(B)
    xlabel('variable');ylabel('coefficient');title('Y-view')
    subplot(2,2,3:4);scatter(U,V)
    xlabel('X-variable');ylabel('Y-variable')
    title('1st pair of canonical variables')
```

Calculate multiple pairs of canonical variables:
```
[A,B,r,U,V] = SCCAec(data1.X,data1.Y,'cxy',[cx,cy],'D',2);

scatter(U(:,1,2),V(:,1,2)) % 2nd pair of canonical variables
```

## PMD-CCA

The `PMDCCA` function works exactly the same way:
```
[A,B,~,~,~,cxy] = PMDCCA(data2.X,data2.Y);
[A,B,r,U,V] = PMDCCA(data1.X,data1.Y,'cxy',[cx,cy]);
[A,B,r,U,V] = PMDCCA(data1.X,data1.Y,'cxy',[cx,cy],'D',2);
```

## StabilityCCA

Run `stabilityCCA` to calculate stability paths and plot them:
```
[A,B] = stabilityCCA(data2.X,data2.Y);

figure
    subplot(2,1,1);plotStabilityCCA(A.probs,A.c,find(data2.A));title('X-view')
    subplot(2,1,2);plotStabilityCCA(B.probs,B.c,find(data2.B));title('Y-view')
```
The ground truth variables are again highlighted.


It is also possible to calculate stability paths for several pairs of canonical variables:
```
[A,B] = stabilityCCA(data1.X,data1.Y,'D',2);

figure
    subplot(2,2,1);plotStabilityCCA(A.probs(:,:,1),A.c,find(data1.A(:,1)))
    title('X-view, 1st pair')
    subplot(2,2,2);plotStabilityCCA(A.probs(:,:,2),A.c,find(data1.A(:,2)))
    title('X-view, 2nd pair')
    subplot(2,2,3);plotStabilityCCA(B.probs(:,:,1),B.c,find(data1.B(:,1)))
    title('Y-view, 1st pair')
    subplot(2,2,4);plotStabilityCCA(B.probs(:,:,2),B.c,find(data1.B(:,2)))
    title('Y-view, 2nd pair')
```

## References
[1] Witten, Daniela M., Robert Tibshirani, and Trevor Hastie. "A penalized matrix decomposition, with applications to sparse principal components and
canonical correlation analysis." Biostatistics 10.3 (2009): 515-534.

[2] Uurtio, Viivi, Sahely Bhadra, and Juho Rousu. "Large-scale sparse kernel canonical correlation analysis." International Conference on Machine Learning.
PMLR, 2019.

[3] Meinshausen, Nicolai, and Peter Bühlmann. "Stability selection." Journal of the Royal Statistical Society: Series B (Statistical Methodology) 72.4
(2010): 417-473.

[4] Shah, Rajen D., and Richard J. Samworth. "Variable selection with error control: another look at stability selection." Journal of the Royal Statistical
Society: Series B (Statistical Methodology) 75.1 (2013): 55-80.
