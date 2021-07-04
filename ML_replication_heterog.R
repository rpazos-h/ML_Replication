
#### Machine Learning Project ####

# This code replicates Fowlie, Holland, & Mansur (2012)
# Then it uses Machine Learning techniques to re-estimate heterogeneous effects

# ML techniques employed:
#     - Honest Trees
#     - GATES & BLP

# Reference
# Fowlie, M., Holland, S. P., & Mansur, E. T. (2012). 
# What do emissions markets deliver and to whom? Evidence 
# from Southern California's NOx trading program. American Economic Review, 102(2), 965-93.

setwd("")

#### Set up ####
# Packages

list.of.packages <-  c('tidyverse', "psych", "stargazer", "data.table",
               "SuperLearner", 'caret',"glmnet", "KernelKnn", "randomForest",
               "reshape2", "xgboost",'kableExtra', 'causalTree','pkgbuild','randomForestCI',
               'causalToolbox', 'forestry', 'grf', 'purrr', 'haven', 'sandwich', 'texreg',
               'lmtest', 'purrr', 'fastDummies', 'haven')

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")

# Load
invisible(lapply(list.of.packages, library, character.only = TRUE))

# Select dplyr
select <- dplyr::select

# Databse

df <- read_dta("https://github.com/rpazos-h/ML_Replication/raw/main/nnmatch_temp.dta")

#### Replication - Matching diff-in-diff ####
##Part to use matched data frame to replicate the matching part of the tables
##First we will try to replicate Table 4 "by hand", this way we could the replicate
##Table 7

df.matched <- read_dta("https://github.com/rpazos-h/ML_Replication/blob/main/match_baseline_3_1.dta?raw=true")

## We filter the data to have correct row ids in the df data

df_na <- df %>% 
  filter(!is.na(DIFFNOX14)& !is.na(PRENOX1))%>% 
  mutate(rowname=1:nrow(.))

## We calculate the weight for the matched data set

weights_ <- df.matched %>% 
  select(id) %>% 
  group_by(id) %>% 
  tally() %>% 
  mutate(weight = 1/ n) %>% 
  select(id, weight)

df.matched  <- df.matched %>%
  left_join(weights_, by = 'id')

# Adjust the treatmente effect for bias due to poor matching

reg_adj <- lm(DIFFNOX14_0 ~ PRENOX1_0m, data = df.matched)
pred_adj0 <- reg_adj$fitted.values
pred_adj1 <- predict(reg_adj, newdata = data.frame(PRENOX1_0m = df.matched$PRENOX1_1m))

# Result
df.matched <- df.matched %>%
  mutate(DIFFNOX14_0_adjust = DIFFNOX14_0 + pred_adj1 - pred_adj0)

##We do a join to get the weight in the dataset

df.matched_filter <- df.matched %>% 
  select(id,index,weight,DIFFNOX14_0_adjust)

df_na_2 <- df_na %>% 
  # Adding weights and Diffnox adjusted
  # Final adjustement of the data
  left_join(df.matched_filter, by = c('rowname' = 'index')) %>% 
  mutate(weight = ifelse(is.na(weight) & dumreclaim==1,1,weight),
         id = ifelse(dumreclaim == 1, rowname, id),
         DIFFNOX14_s = ifelse(dumreclaim==1, DIFFNOX14, DIFFNOX14_0_adjust))

#### TABLE 4 ####
# Levels

reg21 <- lm(DIFFNOX14_s ~ dumreclaim + factor(id), data=df_na_2, weights=weight)

serrors <- list(coeftest(reg21, vcov = vcovHC(reg21, cluster= df_na_2$ab, type = 'HC3'))[,2])
p.values <- list(coeftest(reg21, vcov = vcovHC(reg21, cluster=df_na_2$ab, type = 'HC3'))[,4])

screenreg(reg21, omit.coef = 'id', digits = 2,
          override.se = serrors, override.pvalues = p.values)

# vcovHC HC 3 is the closest to Stata

#### Heterogenous effects ####

##We do the corretion to get the results of the paper
econ_vars_diff <- read_dta('data/econ_vars_diff.dta')
toxic <- read_dta('data/toxic.dta')
#ejmain <- read_dta('data/ej.dta')%>% 
#  select(ab, facid, g_r2009, g_cstl, black1) %>% 
 # .[!duplicated(.),] %>% 
#  mutate(black1= round(black1, digits = 8))

df_na_3 <- df_na_2 %>%
  #mutate(black_ad = round(black1, digits = 6)) %>% 
  #Econ vars diff
  left_join(econ_vars_diff, by = 'co') %>% 
  #Toxic vars
  left_join(toxic, by = c('co', 'ab', 'facid', 'dis')) %>%
  #left_join(ejmain, by = c('facid', 'ab', 'black_ad'= 'black1')) %>% 
  mutate(# Demographics pct of group over total
    #g_r2009 = ifelse(is.na(g_r2009), 0, g_r2009),
    #g_cstl = ifelse(is.na(g_cstl), 0, g_cstl), 
    pctblack1 = 100*(black1)/total1,
    pcthispanic1 = 100*(hispanic1)/total1,
    pctasian1 = 100*(asian1)/total1,
    pctwhite1 = 100*(white1)/total1,
    pctminor1 = 100*(black1 + hispanic1)/total1,
    
    # Pct employment by sector and county
    pctcon = EMP_CON / EMP_ALL,
    pctwhl = EMP_WHL / EMP_ALL,
    pctman = EMP_MAN / EMP_ALL,
    pctret = EMP_RET / EMP_ALL,
    
    # Interactions
    trtpctminor1 = dumreclaim * pctminor1,
    trtPRENOX1 = dumreclaim * PRENOX1,
    trtincome1 = dumreclaim * income1) %>% 
  group_by(dumreclaim) %>% 
  mutate(# Taking Mean by group
         PRENOX1_mean = mean(PRENOX1, na.rm = T),
         income1_mean = mean(income1, na.rm = T),
         pctminor1_mean = mean(pctminor1, na.rm = T),
         pctblack1_mean = mean(pctblack1, na.rm = T),
         pcthispanic1_mean = mean(pcthispanic1, na.rm = T),
         pctasian1_mean = mean(pctasian1, na.rm = T),
         pctwhite1_mean = mean(pctwhite1, na.rm = T),
         trtPRENOX1_mean = mean(trtPRENOX1, na.rm = T),
         trtincome1_mean = mean(trtincome1, na.rm = T),
         trtpctminor1_mean = mean(trtpctminor1, na.rm = T),
         
         # Variable - Mean
         PRENOX1_m = PRENOX1 - PRENOX1_mean,
         income1_m = income1 - income1_mean,
         pctminor1_m = pctminor1 - pctminor1_mean,
         pctblack1_m = pctblack1 - pctblack1_mean,
         pcthispanic1_m = pcthispanic1 - pcthispanic1_mean,
         pctasian1_m = pctasian1 - pctasian1_mean,
         pctwhite1_m = pctwhite1 - pctwhite1_mean,
         trtpctminor1_m = trtpctminor1 - trtpctminor1_mean,
         trtPRENOX1_m = trtPRENOX1 - trtPRENOX1_mean,
         trtincome1_m = trtincome1 - trtincome1_mean) %>% 
  ungroup()


demeans <- function(variables,tr,data){
  dumreclaim <- tr
  names_var <- names(variables)
  df <- data %>% 
    select(all_of(variables), dumreclaim, DIFFNOX14, id, weight, ab) %>% 
    drop_na() %>% 
    group_by(dumreclaim) %>% 
    mutate(across(all_of(variables), ~ .x - mean(.x, na.rm = TRUE), .names = "{.col}_m")) %>% 
    ungroup()
return(df)
}

# Panel A

data11 <- demeans(c('PRENOX1'),'dumreclaim',df_na_3)
reg11_tbl7 <- lm(DIFFNOX14 ~ dumreclaim + PRENOX1 + I(dumreclaim*PRENOX1_m) + factor(id),
             data=data11,weights=weight)

data12 <- demeans(c('PRENOX1', 'income1'),'dumreclaim',df_na_3)
reg12_tbl7 <- lm(DIFFNOX14 ~ dumreclaim + PRENOX1 + I(dumreclaim*income1_m) + income1 + factor(id),
                 data=data12, weights=weight)

data13 <- demeans(c('PRENOX1', 'pctminor1'),'dumreclaim',df_na_3)
reg13_tbl7 <- lm(DIFFNOX14 ~ dumreclaim + PRENOX1 + I(dumreclaim*pctminor1_m) + pctminor1 + factor(id),
                 data=data13,weights=weight)

data14 <- demeans(c('PRENOX1', 'income1'),'dumreclaim',df_na_3)
reg14_tbl7 <- lm(DIFFNOX14 ~ dumreclaim + PRENOX1 + I(dumreclaim*PRENOX1_m) + I(dumreclaim*income1_m) + income1 + 
                   factor(id),
                 data = data14, weights=weight)

data15 <- demeans(c('PRENOX1', 'pctminor1'),'dumreclaim',df_na_3)
reg15_tbl7 <- lm(DIFFNOX14 ~ dumreclaim + PRENOX1 + I(dumreclaim*PRENOX1_m) + I(dumreclaim*pctminor1_m) + pctminor1 +
                       factor(id),
                     data=data15,weights=weight)

data16<-demeans(c('PRENOX1', 'income1', 'pctminor1'),'dumreclaim',df_na_3)
reg16_tbl7 <- lm(DIFFNOX14 ~ dumreclaim + PRENOX1 + I(dumreclaim*income1_m) + I(dumreclaim*pctminor1_m) + pctminor1 +
                   income1 + factor(id),
                 data=data16,weights=weight)

data17<-demeans(c('PRENOX1', 'income1', 'pctminor1'),'dumreclaim',df_na_3)
reg17_tbl7 <- lm(DIFFNOX14 ~ dumreclaim + PRENOX1 + I(dumreclaim*income1_m) + I(dumreclaim*pctminor1_m) + pctminor1 +
                   I(dumreclaim*PRENOX1_m) + income1 + factor(id),
                 data = data17,weights=weight)


serrors <- list(coeftest(reg11_tbl7, vcov = vcovCL(reg11_tbl7, cluster= data11$ab))[,2],
                coeftest(reg12_tbl7, vcov = vcovCL(reg12_tbl7, cluster= data12$ab))[,2],
                coeftest(reg13_tbl7, vcov = vcovCL(reg13_tbl7, cluster= data13$ab))[,2],
                coeftest(reg14_tbl7, vcov = vcovCL(reg14_tbl7, cluster= data14$ab))[,2],
                coeftest(reg15_tbl7, vcov = vcovCL(reg15_tbl7, cluster= data15$ab))[,2],
                coeftest(reg16_tbl7, vcov = vcovCL(reg16_tbl7, cluster= data16$ab))[,2],
                coeftest(reg17_tbl7, vcov = vcovCL(reg17_tbl7, cluster= data17$ab))[,2])

p.values <- list(coeftest(reg11_tbl7, vcov = vcovCL(reg11_tbl7, cluster= data11$ab))[,4],
                 coeftest(reg12_tbl7, vcov = vcovCL(reg12_tbl7, cluster= data12$ab))[,4],
                 coeftest(reg13_tbl7, vcov = vcovCL(reg13_tbl7, cluster= data13$ab))[,4],
                 coeftest(reg14_tbl7, vcov = vcovCL(reg14_tbl7, cluster= data14$ab))[,4],
                 coeftest(reg15_tbl7, vcov = vcovCL(reg15_tbl7, cluster= data15$ab))[,4],
                 coeftest(reg16_tbl7, vcov = vcovCL(reg16_tbl7, cluster= data16$ab))[,4],
                 coeftest(reg17_tbl7, vcov = vcovCL(reg17_tbl7, cluster= data17$ab))[,4])

screenreg(list(reg11_tbl7,reg12_tbl7,reg13_tbl7,reg14_tbl7,reg15_tbl7,reg16_tbl7,reg17_tbl7), digits = 2,
          omit.coef = "(id)|(Intercept)",override.se = serrors, override.pvalues = p.values,
          custom.coef.names = c('Treatment', 'Period 1 NOx', 'Treat x Period 1 NOx', 'Treat x income',
                                'Income', 'Treat x % Minority','% Minority'),
          reorder.coef = c(1,3,4,6,2,5,7),
          custom.model.names = c('(1)','(2)','(3)','(4)','(5)','(6)','(7)'))


texreg(list(reg11_tbl7,reg12_tbl7,reg13_tbl7,reg14_tbl7,reg15_tbl7,reg16_tbl7,reg17_tbl7), digits = 2,
       omit.coef = "(id)|(Intercept)",override.se = serrors, override.pvalues = p.values,
       custom.coef.names = c('Treatment', 'Period 1 NOx', 'Treat x Period 1 NOx', 'Treat x income',
                             'Income', 'Treat x % Minority','% Minority'),
       reorder.coef = c(1,3,4,6,2,5,7),
       custom.model.names = c('(1)','(2)','(3)','(4)','(5)','(6)','(7)'))

#### With machine learning techniques ####

df_na_4 <- df_na_3 %>% 
  select(-rowname:-DIFFNOX14_s) %>% 
  .[!duplicated(.),]

het_df <- df_na_4 %>% 
  ungroup() %>% 
  select(dumreclaim, DIFFNOX14, PRENOX1, PRENOX1_m, income1, income1_m, pctminor1, pctminor1_m,
         pctblack1, pctblack1_m, pcthispanic1, pcthispanic1_m, pctasian1, pctasian1_m, pctwhite1, 
         pctwhite1_m, co, pctcon, pctret, pctwhl, pctman, fsic14, toxic, pdemp, pdpay, pdest) %>% 
  .[complete.cases(.),]

Tr <- het_df[['dumreclaim']]
Y <- het_df[['DIFFNOX14']]

covariates_names <- names(het_df %>% select(PRENOX1, income1,pctblack1, 
                                            pcthispanic1, pctasian1, pdemp,
                                            pdpay, pdest, toxic))

covariates_names_m <- names(het_df %>% select(PRENOX1_m, income1_m,
                                            pctblack1_m, pcthispanic1_m, 
                                            pctasian1_m, pdemp,
                                            pdpay, pdest, toxic))

covariates_paper <- names(het_df %>% select(PRENOX1_m, income1_m, pctminor1_m))

X <- het_df %>% 
  select(all_of(covariates_names), all_of(covariates_names_m), fsic14) %>% 
  # Create one dummy for each variable or put it as factor
  dummy_cols(select_columns = c('fsic14')) 

X_paper <- het_df %>% 
  select(all_of(covariates_paper), fsic14) %>% 
  # Create one dummy for each variable or put it as factor
  dummy_cols(select_columns = c('fsic14')) 

dummies_fsic <-  names(X %>% select(fsic14_1311:fsic14_9711))
#dummies_co <-  names(X %>% select(co_4:co_58))

set.seed(36)
# Split data into 3 samples
folds = createFolds(1:(nrow(het_df)), k=2, list = TRUE) # random split

Y1 <- as.numeric(unlist(Y[folds[[1]]]))
Y2 <- as.numeric(unlist(Y[folds[[2]]]))

Tr1 <- as.numeric(unlist(Tr[folds[[1]]]))
Tr2 <- as.numeric(unlist(Tr[folds[[2]]]))

X1 <- X[folds[[1]],]
X2 <- X[folds[[2]],]

X1_paper <- X_paper[folds[[1]],]
X2_paper <- X_paper[folds[[2]],]

# Creates a vector of 0s and a vector of 1s of length n (hack for later usage)
zeros <- function(n) {
  return(integer(n))
}
ones <- function(n) {
  return(integer(n)+1)
}

###### Honest Trees ######

# Honest Tree 1

honest_tree_formula <- as.formula(paste("Y ~ ", paste(c(covariates_names_m, dummies_fsic), collapse = "+")))

# Training honest causal tree

honest_tree <- honest.causalTree(formula = honest_tree_formula,
                                 data = tibble(Y=Y1,X1[,c(covariates_names_m, dummies_fsic)]),
                                 treatment = Tr1,
                                 est_data = tibble(Y=Y2,X2[,c(covariates_names_m, dummies_fsic)]),
                                 est_treatment = Tr2,
                                 split.Rule = 'CT',
                                 split.Honest = TRUE,
                                 split.Bucket = TRUE,
                                 cv.option = 'CT',
                                 cv.Honest = TRUE,
                                 split.alpha = 0.5,
                                 cv.alpha = 0.5,
                                 bucketMax = 20,
                                 bucketNum = 20,
                                 minsize = 10)

# Plot 
rpart.plot(honest_tree, roundint = F)

## Standard errors
data_1 <- tibble(Y=Y1, X1[,c(covariates_names_m, dummies_fsic)])
data_2 <- tibble(Y=Y2, X2[,c(covariates_names_m, dummies_fsic)])
                                  
######  Standard errors with Honest Tree #######

# On training
honest_tree_pred_1 <- predict(honest_tree, newdata = data_1, 
                              type = "vector")
# On 2  
honest_tree_pred_2 <- predict(honest_tree, newdata = data_2, 
                              type = "vector")

# Adding leaf in samples
data_1$leaf <- as.factor(predict(honest_tree,newdata = data_1,
                                       type = "vector"))

data_2$leaf <- as.factor(predict(honest_tree, newdata = data_2,
                                       type = "vector"))

# lm and standard errors
honest_ols_1 <- lm( Y ~ leaf + Tr * leaf - Tr -1, 
                    data = tibble(Tr = Tr1, data_1))
honest_ols_2 <- lm( Y ~ leaf + Tr * leaf - Tr -1, 
                    data = tibble(Tr = Tr2, data_2))

stargazer(list(honest_ols_1,honest_ols_2), type='text')

## leaf-151.362164109278:Tr Hetereogeneity

##### Sorted Group Average Treatment Effects (GATES) ####
set.seed(1)
gates <- function(Y, W, X, Q=4, prop_scores=F, tree = T) {
  ### STEP 1: split the dataset into two sets, 1 and 2 (50/50)
  split <- createFolds(1:length(Y), k=2)[[1]]
  
  Ya = Y[split]
  Yb = Y[-split]
  
  Xa = X[split]
  Xb = X[-split]
  
  Wa = W[split, ]
  Wb = W[-split, ]
  
  ### STEP 2a: (Propensity score) On set A, train a model to predict X using W. Predict on set B.
  if (prop_scores==T) {
    sl_w1 = SuperLearner(Y = Xa, 
                         X = Wa, 
                         newX = Wb, 
                         family = binomial(), 
                         SL.library = "SL.glm", 
                         cvControl = list(V=0))
    
    p <- sl_w1$SL.predict
  } else {
    p <- rep(mean(Xa), length(Xb))
  }
  
  ### STEP 2b let D = W(set B) - propensity score.
  D <- Xb-p
  
  ### STEP 3a: Get CATE (for example using xgboost) on set A. Predict on set B.
  if(tree == T){
    # Get formula
    tree_fml <- as.formula(paste("Y", paste(names(Wa), collapse = ' + '), sep = " ~ "))
    
    # Causal tree
    honest_tree <- honest.causalTree(formula = tree_fml,
                                     data = data.frame(Y = Ya, Wa),
                                     treatment = Xa,
                                     est_data = data.frame(Y = Yb, Wb),
                                     est_treatment = Xb,
                                     split.alpha = 0.5,
                                     split.Rule = "CT",
                                     split.Honest = T,
                                     cv.alpha = 0.5,
                                     cv.option = "CT",
                                     cv.Honest = T,
                                     split.Bucket = T,
                                     bucketNum = 30,
                                     bucketMax = 10, 
                                     minsize = 30) 
    # Prune
    #opcpid <- which.min(honest_tree$cp[, 4]) 
    #opcp <- honest_tree$cp[opcpid, 1]
    #honest_tree_prune <- prune(honest_tree, cp = opcp)
    
    # CATE
    cate <- predict(honest_tree, newdata = data.frame(Y=Yb, Wb), type = "vector")
    
  } else{
  sl_y = SuperLearner(Y = Ya, 
                      X = data.frame(X=Xa, Wa), 
                      family = gaussian(), 
                      SL.library = "SL.xgboost",
                      cvControl = list(V=0))
  
  pred_0s <- predict(sl_y, data.frame(X=zeros(nrow(Wb)), Wb), onlySL = T)
  pred_1s <- predict(sl_y, newdata=data.frame(X=ones(nrow(Wb)), Wb), onlySL = T)
  
  cate <- pred_1s$pred - pred_0s$pred}
  
  ### STEP 3b: divide the cate estimates into Q tiles, and call this object G. 
  # Divide observations into n tiles
  G <- data.frame(cate) %>% # replace cate with the name of your predictions object
    ntile(Q) %>%  # Divide observations into Q-tiles
    factor()
  
  ### STEP 4: Create a dataframe with Y, W (set B), D, G and p. Regress Y on group membership variables and covariates. 
  df <- data.frame(Y=Yb, Wb, D, G, p)
  
  Wnames <- paste(colnames(Wb), collapse="+")
  fml <- paste("Y ~",Wnames,"+ D:G")
  model <- lm(fml, df, weights = 1/(p*(1-p))) 
  
  clan <- data.frame(Wb, G) %>%
  group_by(G) %>%
  summarise(across(everything(), list(mean = mean, CI_hi = ~ mean(.) + sd(.)*1.645, CI_lo = ~ mean(.) - sd(.)*1.645)))
  
  ### Return list
  result <- list(model = model, clan = clan)
  
  return(result) 
}
  
table_from_gates <-function(model) {
  thetahat <- model$model %>% 
    .$coefficients %>%
    .[c("D:G1","D:G2","D:G3","D:G4")]
  
  # Confidence intervals
  cihat <- confint(model$model)[c("D:G1","D:G2","D:G3","D:G4"),]
  
  res <- tibble(coefficient = c("gamma1","gamma2","gamma3","gamma4"),
                estimates = thetahat,
                ci_lower_90 = cihat[,1],
                ci_upper_90 = cihat[,2])
  
  res_clan <- cbind(res, model$clan)
  
  return(res_clan)
}

output_gates <- rerun(50, table_from_gates(gates(Y = Y, W = X[,c(covariates_names_m, dummies_fsic)], 
                                                 X= Tr, prop_scores = F, tree=F))) %>% # Increase reruns in practice!
  bind_rows %>%
  group_by(coefficient) %>%
  select(-G) %>%
  summarize_all(median)

output_gates[,c(1:4)]

kable(output_gates[,c(1:4)], format='latex', digits = 3,  booktabs = TRUE) %>% 
  kable_minimal()


##### BLP #####

blp <- function(Y, W, X, prop_scores=F) {
  
  ### STEP 1: split the dataset into two sets, 1 and 2 (50/50)
  split <- createFolds(1:length(Y), k=2)[[1]]
  
  Ya = Y[split]
  Yb = Y[-split]
  
  Xa = X[split]
  Xb = X[-split]
  
  Wa = W[split, ]
  Wb = W[-split, ]
  
  ### STEP 2a: (Propensity score) On set A, train a model to predict X using W. Predict on set B.
  if (prop_scores==T) {
    sl_w1 = SuperLearner(Y = Xa, 
                         X = Wa, 
                         newX = Wb, 
                         family = binomial(), 
                         SL.library = "SL.xgboost", 
                         cvControl = list(V=0))
    
    p <- sl_w1$SL.predict
  } else {
    p <- rep(mean(Xa), length(Xb))
  }
  
  ### STEP 2b let D = W(set B) - propensity score.
  D <- Xb-p
  
  ### STEP 3a: Get CATE (for example using xgboost) on set A. Predict on set B.
  sl_y = SuperLearner(Y = Ya, 
                      X = data.frame(X=Xa, Wa), 
                      family = gaussian(), 
                      SL.library = "SL.xgboost", 
                      cvControl = list(V=0))
  
  pred_0s <- predict(sl_y, data.frame(X=zeros(nrow(Wb)), Wb), onlySL = T)
  pred_1s <- predict(sl_y, data.frame(X=ones(nrow(Wb)), Wb), onlySL = T)
  
  cate <- pred_1s$pred - pred_0s$pred
  
  ### STEP 3b: Subtract the expected CATE from the CATE
  C = cate-mean(cate)
  
  ### STEP 4: Create a dataframe with Y, W (set B), D, C and p. Regress Y on W, D and D*C. 
  df <- data.frame(Y=Yb, Wb, D, C, p)
  
  Wnames <- paste(colnames(Wb), collapse="+")
  fml <- paste("Y ~",Wnames,"+ D + D:C")
  model <- lm(fml, df, weights = 1/(p*(1-p))) 
  
  return(model) 
}

table_from_blp <-function(model) {
  thetahat <- model%>% 
    .$coefficients %>%
    .[c("D","D:C")]
  
  # Confidence intervals
  cihat <- confint(model)[c("D","D:C"),]
  
  res <- tibble(coefficient = c("beta1","beta2"),
                estimates = thetahat,
                ci_lower_90 = cihat[,1],
                ci_upper_90 = cihat[,2])
  
  return(res)
}

output_blp <- rerun(50, table_from_blp(blp(Y, X[,c(covariates_names_m)], 
                                           Tr, prop_scores = T))) %>% # Increase reruns in practice!
  bind_rows %>%
  group_by(coefficient) %>%
  summarize_all(median)

output_blp

kable(output_blp, format='latex', digits = 3,  booktabs = TRUE) %>% 
  kable_minimal()

### Exercise replication with matching ####

gates_2 <- function(Y, W, X, weights, Q=4, prop_scores=F, tree = T) {
  ### STEP 1: split the dataset into two sets, 1 and 2 (50/50)
  split <- createFolds(1:length(Y), k=2)[[1]]
  
  Ya = Y[split]
  Yb = Y[-split]
  
  Xa = X[split]
  Xb = X[-split]
  
  Wa = W[split, ]
  Wb = W[-split, ]
  
  weights_a = weights[split]
  weights_b = weights[-split]
  
  ### STEP 2a: (Propensity score) On set A, train a model to predict X using W. Predict on set B.
  if (prop_scores==T) {
    sl_w1 = SuperLearner(Y = Xa, 
                         X = Wa, 
                         newX = Wb, 
                         family = binomial(), 
                         SL.library = "SL.xgboost", 
                         cvControl = list(V=0),
                         obsWeights = weights_a)
    
    p <- sl_w1$SL.predict
  } else {
    p <- rep(mean(Xa), length(Xb))
  }
  
  ### STEP 2b let D = W(set B) - propensity score.
  D <- Xb-p
  
  ### STEP 3a: Get CATE (for example using xgboost) on set A. Predict on set B.
  if(tree == T){
    # Get formula
    tree_fml <- as.formula(paste("Y", paste(names(Wa), collapse = ' + '), sep = " ~ "))
    
    # Causal tree
    honest_tree <- honest.causalTree(formula = tree_fml,
                                     data = data.frame(Y = Ya, Wa),
                                     treatment = Xa,
                                     est_data = data.frame(Y = Yb, Wb),
                                     est_treatment = Xb,
                                     split.alpha = 0.5,
                                     split.Rule = "CT",
                                     split.Honest = T,
                                     cv.alpha = 0.5,
                                     cv.option = "CT",
                                     cv.Honest = T,
                                     split.Bucket = T,
                                     bucketNum = 30,
                                     bucketMax = 10, 
                                     minsize = 40) 
    # Prune
    #opcpid <- which.min(honest_tree$cp[, 4]) 
    #opcp <- honest_tree$cp[opcpid, 1]
    #honest_tree_prune <- prune(honest_tree, cp = opcp)
    
    # CATE
    cate <- predict(honest_tree, newdata = data.frame(Y=Yb, Wb), type = "vector")
    
  } else{
    sl_y = SuperLearner(Y = Ya, 
                        X = data.frame(X=Xa, Wa), 
                        family = gaussian(), 
                        SL.library = "SL.xgboost",
                        cvControl = list(V=0),
                        obsWeights = weights_a)
    
    pred_0s <- predict(sl_y, data.frame(X=zeros(nrow(Wb)), Wb), onlySL = T)
    pred_1s <- predict(sl_y, newdata=data.frame(X=ones(nrow(Wb)), Wb), onlySL = T)
    
    cate <- pred_1s$pred - pred_0s$pred}
  
  ### STEP 3b: divide the cate estimates into Q tiles, and call this object G. 
  # Divide observations into n tiles
  G <- data.frame(cate) %>% # replace cate with the name of your predictions object
    ntile(Q) %>%  # Divide observations into Q-tiles
    factor()
  
  ### STEP 4: Create a dataframe with Y, W (set B), D, G and p. Regress Y on group membership variables and covariates. 
  df <- data.frame(Y=Yb, Wb, D, G, p)
  
  Wnames <- paste(colnames(Wb), collapse="+")
  fml <- paste("Y ~",Wnames,"+ D:G")
  model <- lm(fml, df, weights = 1/(p*(1-p))) 
  
  clan <- data.frame(Wb, G) %>%
    group_by(G) %>%
    summarise(across(everything(), list(mean = mean, CI_hi = ~ mean(.) + sd(.)*1.645, CI_lo = ~ mean(.) - sd(.)*1.645)))
  
  ### Return list
  result <- list(model = model, clan = clan)
  
  return(result) 
}

het_wmatching <- tibble('Y' = df_na_3$DIFFNOX14, 'Tr' = df_na_3$dumreclaim, 
                        df_na_3[,c(covariates_paper, 'id', 'fsic14')], 
                        'weights' = df_na_3$weight) %>% 
  drop_na(.) %>% 
  dummy_cols(select_columns = c("fsic14", "id"))

dummies_id <- names(het_wmatching %>% select(id_318:id_1069))

output_gates_matching <- rerun(10,table_from_gates(gates_2(Y = het_wmatching$Y, W = het_wmatching[,c(covariates_paper,dummies_id)], 
                                                           X= het_wmatching$Tr, prop_scores = F, tree = F, weights = het_wmatching$weights))) %>% 
  bind_rows %>%
  group_by(coefficient) %>%
  select(-G) %>%
  summarize_all(median)

output_gates_matching[,c(1:4)]

##### BLP #####

blp_2 <- function(Y, W, X, weights, prop_scores=F) {
  
  ### STEP 1: split the dataset into two sets, 1 and 2 (50/50)
  split <- createFolds(1:length(Y), k=2)[[1]]
  
  Ya = Y[split]
  Yb = Y[-split]
  
  Xa = X[split]
  Xb = X[-split]
  
  Wa = W[split, ]
  Wb = W[-split, ]
  
  weights_a = weights[split]
  weights_b = weights[-split]
  
  ### STEP 2a: (Propensity score) On set A, train a model to predict X using W. Predict on set B.
  if (prop_scores==T) {
    sl_w1 = SuperLearner(Y = Xa, 
                         X = Wa, 
                         newX = Wb, 
                         family = binomial(), 
                         SL.library = "SL.xgboost", 
                         cvControl = list(V=0),
                         obsWeights = weights_a)
    
    p <- sl_w1$SL.predict
  } else {
    p <- rep(mean(Xa), length(Xb))
  }
  
  ### STEP 2b let D = W(set B) - propensity score.
  D <- Xb-p
  
  ### STEP 3a: Get CATE (for example using xgboost) on set A. Predict on set B.
  sl_y = SuperLearner(Y = Ya, 
                      X = data.frame(X=Xa, Wa), 
                      family = gaussian(), 
                      SL.library = "SL.xgboost", 
                      cvControl = list(V=0),
                      obsWeights = weights_a)
  
  pred_0s <- predict(sl_y, data.frame(X=zeros(nrow(Wb)), Wb), onlySL = T)
  pred_1s <- predict(sl_y, data.frame(X=ones(nrow(Wb)), Wb), onlySL = T)
  
  cate <- pred_1s$pred - pred_0s$pred
  
  ### STEP 3b: Subtract the expected CATE from the CATE
  C = cate-mean(cate)
  
  ### STEP 4: Create a dataframe with Y, W (set B), D, C and p. Regress Y on W, D and D*C. 
  df <- data.frame(Y=Yb, Wb, D, C, p)
  
  Wnames <- paste(colnames(Wb), collapse="+")
  fml <- paste("Y ~",Wnames,"+ D + D:C")
  model <- lm(fml, df, weights = 1/(p*(1-p))) 
  
  return(model) 
}

output_blp <- rerun(50, table_from_blp(blp_2(Y = het_wmatching$Y, W = het_wmatching[,c(covariates_paper,dummies_id)], 
                                             X= het_wmatching$Tr, prop_scores = F, weights = het_wmatching$weights))) %>% # Increase reruns in practice!
  bind_rows %>%
  group_by(coefficient) %>%
  summarize_all(median)

output_blp

kable(output_blp, format='latex', digits = 3,  booktabs = TRUE) %>% 
  kable_minimal()


