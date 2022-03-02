#install.packages('survivalROC')
library(Hmisc)
library(PredictABEL)  #
library(foreign)
library(survival)
library(rms)  #
library(dplyr)
library(survival)
library(pec)
library(rpart) #classification and regression trees
library(partykit) #treeplots
library(MASS) #breast and pima indian data
library(ElemStatLearn) #prostate data
library(randomForest) #random forests
library(xgboost) #gradient boosting 
library(caret) #tune hyper-parameters #
library(regplot)
library(mstate)
library(survivalROC)
library(ggplot2)
library("DynNom", lib.loc="~/R/win-library/3.6") #
#install.packages('rsconnect')    
library(rsconnect)
library(DynNom)    #
library("prodlim")
library("survcomp")

#install.packages('tidyselect')


setwd("E:\\")

#GI_bleed<-read.csv("GI_bleed.csv") 

##将数据集分为训练集跟测试集2:8
#train_1 <-createDataPartition(y=GI_bleed$patientunitstayid,p=0.8,list=FALSE)
#train <- GI_bleed[train_1, ]
#test <- GI_bleed[-train_1, ]
#write.csv(train, 'train.csv')
#write.csv(test, 'test.csv')

#train1 <- data.frame(reduce_redundency(train, threshold = 0.9)$dat.redd) #去冗余
#write.csv(train1, "train1.csv")
##生存资料lasso回归筛选变量
library(glmnet)
library(survival)
library(foreign)
GI_bleed<-read.csv("train.csv")
GI_bleed<-as.data.frame(GI_bleed)
GI_bleed<-na.omit(GI_bleed)
head(GI_bleed)
x <- as.matrix(GI_bleed[,9:59])
y <- Surv(GI_bleed$icu_los_day,GI_bleed$live_status==1)
lasso <- glmnet(x, y, family = "cox", alpha = 1)
print(lasso)
plot(lasso, xvar = "lambda", label = TRUE)

##交叉验证Lasso回归
#type.measure=deviance 使用deviance，即-2倍的Log-likelihood
#type.measure=mse 使用拟合因变量与实际应变量的mean squred error
#type.measure=mae 使用mean absolute error
#type.measure=class 使用模型分类的错误率(missclassification error)
#type.measure=auc or C (FOR SURVIVAL) 使用area under the ROC curve，是现在最流行的综合考量模型性能的一种参数

set.seed(123)
fitCV <- cv.glmnet(x, y, family = "cox",
                   type.measure = "deviance",
                   nfolds = 50)

plot(fitCV)

fitCV$lambda.1se
coef(fitCV, s = "lambda.1se")

#fitCV$lambda.min
#coef(fitCV, s = "lambda.min")

#fitCV$lambda.min
#coef(fitCV, s = "lambda.min")

##建模
train<-read.csv("train.csv") 
test<-read.csv("test.csv") 
test1<-read.csv("test1.csv") 
test2<-read.csv("test2.csv") 
test_sepsis<-read.csv("eicu_sepsis.csv") 
head(train)
dd=datadist(train)
option <- options(datadist = "ddist")
options(datadist="dd") 
attach(train)
as.matrix(train)
str(train)


attach(test)
str(test)
as.matrix(test)

attach(test1)
str(test1)
as.matrix(test1)

attach(test2)
str(test2)
as.matrix(test2)

attach(test_sepsis)
str(test_sepsis)
as.matrix(test_sepsis)
##
Srv = Surv(train$icu_los_day, train$live_status)


##30天死亡lasso-cox
as.double(inr_max)
as.numeric(albumin_min)
as.numeric(potassium )
as.integer(oobventday1)
as.integer(vasopressor)

coxm1 <- cph(Surv(icu_los_day,live_status)~age_65+oobventday1+rr+Spo2+sbp+inr_max
             +bicarbonate_min+albumin_min+potassium+vasopressor, data=train,surv=T)
#coxm1
#summary(coxm1)
DynNom(coxm1) ##模型可视化
DNbuilder(coxm1)

surv <- Survival(coxm1)
surv1 <- function(x) (1-surv(14, x))
surv2 <- function(x) (1-surv(28, x))

nom <- nomogram(coxm1, fun=list(surv1,surv2),lp = F,
                funlabel=c("14-day Death Probability", "28-day Death Probability"),
                maxscale=100,
                fun.at=c('0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9'))
plot(nom, xfrac=.2)

nom ###计算分数



##使用方差膨胀因子(VIF)检验多重共线性。
#当0<VIF<10，不存在多重共线性（注意：在《R语言实战》
#第2版P182中认为VIF>4就存在多重共线性）；
#当10≤VIF<100，存在较强的多重共线性，当VIF>=100，多重共线性非常严重
library(car)
vif(coxm1)

##计算points
as.matrix(test)
as.matrix(test1)
library(nomogramFormula)
options(option)
results <- formula_rd(nomogram=nom)
train$points <- points_cal(formula = results$formula,rd=train)
test$points <- points_cal(formula = results$formula,rd=test)
test1$points <- points_cal(formula = results$formula,rd=test1)
test2$points <- points_cal(formula = results$formula,rd=test2)
test_sepsis$points <- points_cal(formula = results$formula,rd=test_sepsis)

write.csv(train,'train.csv')
write.csv(test,'test.csv')
write.csv(test1,'test1.csv')
write.csv(test2,'test2.csv')
write.csv(test_sepsis,'test_sepsis.csv')


###计算预测概率
##train
t <- c(14)
survprob14 <- predictSurvProb(coxm1,newdata=train,times=t)
head(survprob14)
train_pre <- cbind(train, survprob14)
write.csv(train_pre, 'train_pre.csv')

t <- c(28)
survprob28 <- predictSurvProb(coxm1,newdata=train_pre,times=t)
head(survprob28)
train_pre <- cbind(train_pre, survprob28)
str(train_pre)
write.csv(train_pre, 'train_pre.csv')

##e-ICU internal test
t <- c(14)
survprob14 <- predictSurvProb(coxm1,newdata=test,times=t)
head(survprob14)
test_pre <- cbind(test, survprob14)
write.csv(test_pre, 'test_pre.csv')

t <- c(28)
survprob28 <- predictSurvProb(coxm1,newdata=test_pre,times=t)
head(survprob28)
test_pre <- cbind(test_pre, survprob28)
str(test_pre)
write.csv(test_pre, 'test_pre.csv')


##test1 Tongji cohort
t <- c(14)
survprob14 <- predictSurvProb(coxm1,newdata=test1,times=t)
head(survprob14)
test1_pre <- cbind(test1, survprob14)
write.csv(test1_pre, 'test1_pre.csv')

t <- c(28)
survprob28 <- predictSurvProb(coxm1,newdata=test1_pre,times=t)
head(survprob28)
test1_pre <- cbind(test1_pre, survprob28)
str(test1_pre)
write.csv(test1_pre, 'test1_pre.csv')

##test2 MIMIC cohort
t <- c(14)
survprob14 <- predictSurvProb(coxm1,newdata=test2,times=t)
head(survprob14)
test2_pre <- cbind(test2, survprob14)
write.csv(test2_pre, 'test2_pre.csv')

t <- c(28)
survprob28 <- predictSurvProb(coxm1,newdata=test2_pre,times=t)
head(survprob28)
test2_pre <- cbind(test2_pre, survprob28)
str(test2_pre)
write.csv(test2_pre, 'test2_pre.csv')

#训练集C-index
library("survival")
library("prodlim")
library("survcomp")
f1<-coxph(Surv(icu_los_day,live_status==1)~age_65+oobventday1+rr+Spo2+sbp+inr_max
          +bicarbonate_min+albumin_min+potassium+vasopressor,data=train)
summary(f1)
sum.surv<-summary(f1)
as.matrix(train)
c_index<-sum.surv$concordance
c_index


f1_1<-coxph(Surv(icu_los_day,live_status==1)~oakland_score,data=train)
summary(f1_1)
sum.surv<-summary(f1_1)
as.matrix(train)
c_index_1<-sum.surv$concordance
c_index_1

#比较方法1
anova(f1,f1_1) 
# 方法2
library(rms)
library(lmtest)
compare <- lrtest(f1, f1_1)
compare

#内部测试 c-index

library(survival)
f2<-coxph(Surv(icu_los_day,live_status==1)~points,data=test)
summary(f2)
sum.surv<-summary(f2)
c_index_1<-sum.surv$concordance
c_index_1

f2_1<-coxph(Surv(icu_los_day,live_status==1)~sofatotal,data=test)
compare <- lrtest(f2, f2_1)
compare

#外部测试 c-index 1 Tongji cohort

library(survival)
f3<-coxph(Surv(icu_los_day,live_status==1)~points,data=test1)
summary(f3)
sum.surv<-summary(f3)
c_index_2<-sum.surv$concordance
c_index_2

f3_1<-coxph(Surv(icu_los_day,live_status==1)~sofatotal,data=test1)
compare <- lrtest(f3, f3_1)
compare


#外部测试 c-index 2 MIMIC

library(survival)
f4<-coxph(Surv(icu_los_day,live_status==1)~points,data=test2)
summary(f4)
sum.surv<-summary(f4)
c_index_3<-sum.surv$concordance
c_index_3


f4_1<-coxph(Surv(icu_los_day,live_status==1)~sofa,data=test2)
compare <- lrtest(f4, f4_1)
compare


#外部测试 c-index sepsis

library(survival)
f5<-coxph(Surv(icu_los_day,live_status==1)~points,data=test_sepsis)
summary(f5)
sum.surv<-summary(f5)
c_index_4<-sum.surv$concordance
c_index_4


anova(f1,f5) 


##训练集校正曲线###################
##14day
s<-Surv(train$icu_los_day,train$live_status,type="right")
f <- cph(s~age_65+oobventday1+rr+Spo2+sbp+inr_max
         +bicarbonate_min+albumin_min+potassium+vasopressor, x=TRUE, y=TRUE,surv = TRUE,time.inc=14,data=train)

cal<-calibrate(f,u=14,cmethod='KM',m=100, method='boot', B=100)

plot(cal,xlim = c(0,1),ylim= c(0,1),
     errbar.col=c(rgb(0,0,0,maxColorValue=255)),col=c(rgb(255,0,0,maxColorValue=255)))

abline(0,1,lty=3,lwd=2,col=c(rgb(0,0,255,maxColorValue= 255)))

summary(cal)

##28day
s<-Surv(train$icu_los_day,train$live_status,type="right")
f <- cph(s~age_65+oobventday1+rr+Spo2+sbp+inr_max
         +bicarbonate_min+albumin_min+potassium+vasopressor, x=TRUE, y=TRUE,surv = TRUE,time.inc=28,data=train)

cal<-calibrate(f,u=28,cmethod='KM',m=100, method='boot', B=100)

plot(cal,xlim = c(0,1),ylim= c(0,1),
     errbar.col=c(rgb(0,0,0,maxColorValue=255)),col=c(rgb(255,0,0,maxColorValue=255)))

abline(0,1,lty=3,lwd=2,col=c(rgb(0,0,255,maxColorValue= 255)))

summary(cal)



#内部验证校正曲线###################
setwd("GI_bleed.csv")
test<-read.csv("test.csv") 
head(test)
dd=datadist(test)
options(datadist="dd") 

##14day
s<-Surv(test$icu_los_day,test$live_status,type="right")
f <- cph(s~points, x=TRUE, y=TRUE,surv = TRUE,time.inc=14,data=test)

cal<-calibrate(f,u=14,cmethod='KM',m=100, method='boot', B=100)

plot(cal,xlim = c(0,1),ylim= c(0,1),
     errbar.col=c(rgb(0,0,0,maxColorValue=255)),col=c(rgb(255,0,0,maxColorValue=255)))

abline(0,1,lty=3,lwd=2,col=c(rgb(0,0,255,maxColorValue= 255)))

##28day
s<-Surv(test$icu_los_day,test$live_status,type="right")
f <- cph(s~points, x=TRUE, y=TRUE,surv = TRUE,time.inc=28,data=test)

cal<-calibrate(f,u=10,cmethod='KM',m=100, method='boot', B=100)

plot(cal,xlim = c(0,1),ylim= c(0,1),
     errbar.col=c(rgb(0,0,0,maxColorValue=255)),col=c(rgb(255,0,0,maxColorValue=255)))

abline(0,1,lty=3,lwd=2,col=c(rgb(0,0,255,maxColorValue= 255)))


##外部校正曲线1 tingji###################
##14day
setwd("GI_bleed.csv")
test1<-read.csv("test1.csv") 
head(test1)
dd=datadist(test1)
options(datadist="dd") 

##14day
s<-Surv(test1$icu_los_day,test1$live_status,type="right")
f <- cph(s~points, x=TRUE, y=TRUE,surv = TRUE,time.inc=14,data=test1)

cal<-calibrate(f,u=10,cmethod='KM',m=100, method='boot', B=100)

plot(cal,xlim = c(0,1),ylim= c(0,1),
     errbar.col=c(rgb(0,0,0,maxColorValue=255)),col=c(rgb(255,0,0,maxColorValue=255)))

abline(0,1,lty=3,lwd=2,col=c(rgb(0,0,255,maxColorValue= 255)))

##28day
s<-Surv(test1$icu_los_day,test1$live_status,type="right")
f <- cph(s~points, x=TRUE, y=TRUE,surv = TRUE,time.inc=28,data=test1)

cal<-calibrate(f,u=10,cmethod='KM',m=100, method='boot', B=100)

plot(cal,xlim = c(0,1),ylim= c(0,1),
     errbar.col=c(rgb(0,0,0,maxColorValue=255)),col=c(rgb(255,0,0,maxColorValue=255)))

abline(0,1,lty=3,lwd=2,col=c(rgb(0,0,255,maxColorValue= 255)))


#外部验证校正曲线2###################
setwd("GI_bleed.csv")
test2<-read.csv("test2.csv") 
head(test2)
dd=datadist(test2)
options(datadist="dd") 

##14day
s<-Surv(test2$icu_los_day,test2$live_status,type="right")
f <- cph(s~points, x=TRUE, y=TRUE,surv = TRUE,time.inc=14,data=test2)

cal<-calibrate(f,u=14,cmethod='KM',m=100, method='boot', B=100)

plot(cal,xlim = c(0,1),ylim= c(0,1),
     errbar.col=c(rgb(0,0,0,maxColorValue=255)),col=c(rgb(255,0,0,maxColorValue=255)))

abline(0,1,lty=3,lwd=2,col=c(rgb(0,0,255,maxColorValue= 255)))

##28day
s<-Surv(test2$icu_los_day,test2$live_status,type="right")
f <- cph(s~points, x=TRUE, y=TRUE,surv = TRUE,time.inc=28,data=test2)

cal<-calibrate(f,u=20,cmethod='KM',m=100, method='boot', B=100)

plot(cal,xlim = c(0,1),ylim= c(0,1),
     errbar.col=c(rgb(0,0,0,maxColorValue=255)),col=c(rgb(255,0,0,maxColorValue=255)))

abline(0,1,lty=3,lwd=2,col=c(rgb(0,0,255,maxColorValue= 255)))


#sepsis校正曲线###################
setwd("D:/GI bleed")
test_sepsis<-read.csv("test_sepsis.csv") 
head(test_sepsis)
dd=datadist(test_sepsis)
options(datadist="dd") 

##14day
s<-Surv(test_sepsis$icu_los_day,test_sepsis$live_status,type="right")
f <- cph(s~points, x=TRUE, y=TRUE,surv = TRUE,time.inc=14,data=test_sepsis)

cal<-calibrate(f,u=14,cmethod='KM',m=100, method='boot', B=100)

plot(cal,xlim = c(0,1),ylim= c(0,1),
     errbar.col=c(rgb(0,0,0,maxColorValue=255)),col=c(rgb(255,0,0,maxColorValue=255)))

abline(0,1,lty=3,lwd=2,col=c(rgb(0,0,255,maxColorValue= 255)))

##28day
s<-Surv(test_sepsis$icu_los_day,test_sepsis$live_status,type="right")
f <- cph(s~points, x=TRUE, y=TRUE,surv = TRUE,time.inc=28,data=test2)

cal<-calibrate(f,u=20,cmethod='KM',m=100, method='boot', B=100)

plot(cal,xlim = c(0,1),ylim= c(0,1),
     errbar.col=c(rgb(0,0,0,maxColorValue=255)),col=c(rgb(255,0,0,maxColorValue=255)))

abline(0,1,lty=3,lwd=2,col=c(rgb(0,0,255,maxColorValue= 255)))


##决策曲线
library(rmda) 

setwd("D:/GI bleed")
GI_bleed<-read.csv("GI_bleed.csv") 
total<-read.csv("total.csv") 
head(total)
attach(total)
str(total)



#决策曲线新代码
#颜色代码"black"   "red"     "green3"  "blue"    "cyan"    "magenta" "yellow"  "gray" 
#模型1 
#cph(Surv(icu_los_day,live_status==1)~oobventday1+rr+Spo2+sbp+inr_max
#    +bicarbonate_min+albumin_min+potassium+vasopressor, data=train,surv=T)

set.seed(123)
nomogram <- lrm(live_status~points,data=total,x=T,y=T)
nomogram <- decision_curve(live_status~points,data=total,
                           
                           thresholds = seq(0, 1.0, by = .005),
                           
                           confidence.intervals =FALSE,study.design = 'case-control',population.prevalence = 0.1)

plot_decision_curve (nomogram, curve.names = "nomogram")

##
sofatotal <- decision_curve(live_status~sofatotal,data=total,
                            
                            thresholds = seq(0, 1.0, by = .005),
                            
                            confidence.intervals =FALSE,study.design = 'case-control',population.prevalence = 0.1)

plot_decision_curve (sofatotal, curve.names = "sofatotal")

##
gbs_score <- decision_curve(live_status~gbs_score,data=total,
                            
                            thresholds = seq(0, 1.0, by = .005),
                            
                            confidence.intervals =FALSE,study.design = 'case-control',population.prevalence = 0.1)

plot_decision_curve (gbs_score, curve.names = "gbs_score")

##
aims65 <- decision_curve(live_status~aims65,data=total,
                            
                            thresholds = seq(0, 1.0, by = .005),
                            
                            confidence.intervals =FALSE,study.design = 'case-control',population.prevalence = 0.1)

plot_decision_curve (aims65, curve.names = "aims65")



##
oakland_score <- decision_curve(live_status~oakland_score,data=total,
                         
                         thresholds = seq(0, 1.0, by = .005),
                         
                         confidence.intervals =FALSE,study.design = 'case-control',population.prevalence = 0.1)

plot_decision_curve (oakland_score, curve.names = "oakland_score")

##
sengupta_score <- decision_curve(live_status~sengupta_score,data=total,
                                
                                thresholds = seq(0, 1.0, by = .005),
                                
                                confidence.intervals =FALSE,study.design = 'case-control',population.prevalence = 0.1)

plot_decision_curve (sengupta_score, curve.names = "sengupta_score")




List<- list(nomogram,sofatotal, gbs_score, aims65,oakland_score,sengupta_score)

plot_decision_curve(List,curve.names= c('nomogram','sofatotal', 'gbs_score', 'aims65',
                                        'oakland_score','sengupta_score'),cost.benefit.axis =FALSE,col = c('red','blue','green3',
                                                                                                            'black','grey','green'),confidence.intervals =FALSE,standardize = FALSE)


#"black"   "red"     "green3"  




set.seed(123)
Nomogram <- decision_curve(status~age + sex + Respiratory_system_disease+ AST_1 + Hs_CRP_1
                           +Hypersensitive_troponin_1	+Leu_1 +Lym_2	+nitrogen_1	+PCT_2,
                           
                           data = test,
                           
                           thresholds = seq(0, 1.0, by = .005),
                           
                           confidence.intervals =FALSE,study.design = 'case-control',population.prevalence = 0.3)

plot_decision_curve (Nomogram, curve.names = "radiomics model")
##age+sex+Hypertension+ Diabetes+Respiratory_system_disease+ Cardiovascular_disease+smoking+Leu_1+Lym_1+Neu_1+PLT_1+ Hb_1+PCT_1+Hypersensitive_troponin_1+Hs_CRP_1+
 ## LDH_1+Serum_creatinine_1+nitrogen_1+TB_1+AST_1+ALT_1+Albumin_1+D_dimer_1+PT_1


##
