
# 1. Cell Forest
## Description
Cell Forest is a method for single cell clustering that is aimed at improving clustering accuracies through feature selection. Based on a large-scale analysis of 20 real world scRNA-seq datasets that cover a wide spectrum of biological scenarios, it was found that Cell Forest achieved, on average, superior performances to six benchmark methods inlcuding RAFSIL2, SIMLR and SC3. 

## Download
Cell Forest is implemented as an R package, which is freely available for non-commercial use. 

Version Changes 
[CellForest_0.2.0.tar.gz](https://github.com/BioinformaticsCSU/CellForest/blob/master/CellForest_0.2.0.tar.gz)

# 2. Install

- Step 1: Firstly, install the dependent packages: **gplots, foreach, doParallel and randomForest**.

- Step 2: Download the above CellForest package and install it in R (tested on version 3.2.0)




# 3. Usage
Notes: CellForest was tested on linux, Mac and Windows; and it runs smoothly on these different systems.

Using CellForest is very simple. Just follow the steps below: 

Step 1: open your R or Rstudio 
Step 2: in the R command window, run the following command to load the R package
```
> library(CellForest)
```
Step 3: in R command window, run the following command to see the help document for running Cell Forest. Then, you should be able to see a help page.
```
> ?CellForest
```
Step 4: At the end of the help page, there is an example code. Copy these codes to command to run as follows:
```
data(CFDemo)
result = CellForest(data,kcluster = kprior)
performance = evalcluster(label,result$cluster)
```

# 4. Contact
If any questions, please do not hesitate to contact us at: 

Hongdong Li, hongdong@csu.edu.cn

Jianxin Wang, jxwang@csu.edu.cn


# 5. How to cite?
If you use this tool, please cite the following work.

Hongdong Li, Yunpei Xu, Cuixiang Lin, Fangxiang Wu and Jianxin Wang, Cell Forest: Accurate Clustering of Single Cells Through Feature Selection, 2018, submitted  
