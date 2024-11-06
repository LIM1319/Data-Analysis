##Q2(a)
## Simulate the data ##
n <- 500
X1 <- runif(n, -10, 10)
X2 <- rnorm(n, mean = 0, sd = 2)
X3 <- rbinom(n, size = 1, prob = 0.7)
Xmat <- cbind(rep(1,n), X1, X2, X3)
b0 <- -0.8
b1 <- 0.1
b2 <- 0.2
b3 <- 0.3
eta0<-(b0+b1*X1+b2*X2+b3*X3)
mu0 <- exp(eta0)
Y <- rpois(n, lambda = exp(eta0))

##Q2(b)
##CREATE A FUNCTION FOR IRLS ##
glm.utsg <- function(Y, Xmat, tol = 1e-8){
  #1 Initialization
  beta <-rep(0, ncol(Xmat))
  eta <- Xmat %*% beta
  mu <- exp(eta)
  #2.Set W^(1)  W=diag(μ1,μ2,...,μn)
  W <- diag(as.vector(mu))
  #3 For i = 1, 2,.., n, set zi^(1)
  Z <- eta + solve(W) %*% (Y - mu)
  #Set ∂ℓ(β)/∂β| β = ˆβ
  eqns <- sum(t(Xmat) %*% (Y-mu))
  istep <- 0
  
  #7 Repeat 2-6 until ∂ℓ(β)/∂β| β = ˆβ < δ
  while(eqns^2 > tol){
    #4 Estimate ˆβ(1)
    beta <- solve(t(Xmat) %*% W %*% Xmat) %*% t(Xmat) %*% W %*% Z
    #5 Set ˆηi(1) and ˆµi(1)
    eta <- Xmat %*% beta
    mu <- exp(eta)
    #6 Perform 2-4 again to obtain W^(2), zi(2) and ˆβ^(2).
    W <- diag(as.vector(mu))
    Z <- eta + solve(W) %*% (Y - mu)
    #Calculate ∂ℓ(β)/∂β for new ˆβ
    eqns <- sum(t(Xmat) %*% (Y - mu))
    istep <- istep +1
  }
  SE <- sqrt(diag(solve(t(Xmat) %*% W %*% Xmat)))
  z = beta/SE
  res.mat <- data.frame(Estimates = beta, OR = c(1, exp(beta)[-1]),
                        SE = SE,
                        z = z,
                        p_value = ifelse(z < 0, 2*pnorm(z), 2*(1 - pnorm(z)))
  )
  rownames(res.mat)[1] <- "Intercept"
  results <- list(Table = res.mat, Iteration = istep)
  return(results)
}  
## run the function and compare it to R's glm() ##
glm.utsg(Y = Y, Xmat = Xmat)
## Run the glm() code in R ##
glm.p <- glm(Y ~ X1 + X2 + X3, family = poisson)
summary(glm.p)

##Q2(c)
## Function to simulate data and check coverage for β3 ##
simulation <- function(n, num_simulations = 10000, beta3 = 0.3) {
  # Initialize counter for successful coverage
  count <- 0
  
  # Loop over the number of simulations
  for (i in 1:num_simulations) {
    
    # Simulate data
    X1 <- runif(n, -10, 10)
    X2 <- rnorm(n, mean = 0, sd = 2)
    X3 <- rbinom(n, size = 1, prob = 0.7)
    Xmat <- cbind(rep(1,n), X1, X2, X3)
    b0 <- -0.8
    b1 <- 0.1
    b2 <- 0.2
    b3 <- 0.3
    eta0<-(b0+b1*X1+b2*X2+b3*X3)
    mu0 <- exp(eta0)
    Y <- rpois(n, lambda = exp(eta0))
    
    # Fit the model 
    model <- glm(Y ~ X1 + X2 + X3, family = poisson(link = "log"))
    
    # Get the Wald confidence interval for β3
    conf_int <- confint.default(model)[4, ]
    
    # Check if the true value of β3 is within the confidence interval
    if (conf_int[1] <= beta3 && conf_int[2] >= beta3) {
      count <- count + 1
    }
  }
  
  # Calculate and return the coverage probability
  probability <- count / num_simulations
  return(probability)
}

# Sample size n = 500
paste("Coverage Probability for n = 500: ", simulation(n = 500))

# Sample size n = 100
paste("Coverage Probability for n = 100: ", simulation(n = 100))

# Sample size n = 30
paste("Coverage Probability for n = 30: ", simulation(n = 30))



#Q3(a)
library(tidyverse)
library(NHANES)
small.nhanes <- NHANES %>% filter(SurveyYr == "2011_12", Age > 17) %>%
  select(1, 3, 4, 8:11, 13, 25, 61) %>% na.omit() %>% group_by(ID) %>%
  slice(1) %>% ungroup()
# Set seed using student ID 
set.seed(1008601452)
# Randomly select 500 observations
sample_data <- small.nhanes %>% sample_n(500)
# Fit a logistic regression model
full_model <- glm(SmokeNow ~ . -ID, data = sample_data, family = binomial)
summary(full_model)


#Q3(b)
# using AIC
stepwise_aic <- step(full_model, trace = 0)
summary(stepwise_aic)

# using BIC
stepwise_bic <- step(full_model, trace = 0, k = log(nrow(sample_data)))
summary(stepwise_bic)

#LASSO
library(glmnet)
X <- model.matrix(full_model)[, -1]
y <- as.numeric(sample_data$SmokeNow)
lasso_model <- cv.glmnet(X, y, family = "binomial" , type.measure = "class")
plot(lasso_model)
coef(lasso_model, s = "lambda.min")
#Determined zero effect.We refit the model excluding these variables.
lasso<- glm(SmokeNow ~  Gender + Age + Poverty + BPSysAve, data = sample_data, family = binomial)
summary(lasso)

anova(stepwise_aic, stepwise_bic)
anova(stepwise_aic, lasso)


#Q3(c)
#model validation
library(rms)
lrm.final <- lrm(SmokeNow ~  Gender + Age+ Poverty + BPSysAve, data=sample_data,
                 x =TRUE, y = TRUE, model= T)
cross.calib <- calibrate(lrm.final, method="crossvalidation", B=10) # model calibration
plot(cross.calib, las=1, xlab = "Predicted Probability")

#Q3(d)
library(pROC)
roc.logit <- roc(lasso$y ~ fitted(lasso))
plot(1-roc.logit$specificities, roc.logit$sensitivities, xlim = c(0,1),
     ylim = c(0,1), type = 'l', lwd = 2, col = 'red', bty = 'n',
     xlab = "1-Specificity", ylab = "Sensitivity")
abline(0, 1, lty = 2, col = 'blue')
text(0.7,0.4,label = paste("AUC = ", round(auc(roc.logit),2)))

#Q3(e)
# Find the remaining 310 observations that were not included in the original sample
remaining_data <- small.nhanes %>% anti_join(sample_data, by = "ID")
lasso2<- glm(SmokeNow ~  Gender + Age + Poverty + BPSysAve, data = remaining_data, family = binomial)
remaining_data <- remaining_data %>% mutate(predprob = fitted(lasso2),linpred = predict(lasso2))

(hldf <- remaining_data %>% group_by(ntile(linpred, 10)) %>%
    summarise(y = sum(SmokeNow == "Yes"), ppred = mean(predprob), count = n()) %>%
    mutate(se.fit = sqrt(ppred*(1-ppred)/count)))

hldf %>%
  ggplot(aes(x = ppred, y = y / count, ymin = y / count - 2 * se.fit, ymax = y / count + 2 * se.fit)) +
  geom_point() +
  geom_linerange(color = grey(0.75)) +
  geom_abline(intercept = 0, slope = 1) 

# computing HL test statistic in R
with(hldf, sum((y-count*ppred)^2/(count*ppred*(1-ppred))))
# critical value
qchisq(0.95, 8)
# p-value
1-pchisq(8.494, 8)



#Q4(a)
library(faraway)
data(hsb)
# Gender and Program Choice
gender_program_table <- table(hsb$gender, hsb$prog)
gender_program_prop <- prop.table(gender_program_table, margin = 1)
gender_program_prop
# SES and Program Choice
ses_program_table <- table(hsb$ses, hsb$prog)
ses_program_prop <- prop.table(ses_program_table, margin = 1)
ses_program_prop


#Q4(b)
library(nnet)
# Fit the multinomial logistic regression model
model <- multinom(prog ~ gender + race + ses + schtyp + read + write + math + science + socst, data = hsb)
summary(model)
# Calculate the odds ratios for math and science
coefficients_vocation <- summary(model)$coefficients["vocation", c("math", "science")]
odds_ratios <- exp(coefficients_vocation)
odds_ratios
# Calculate the 95% Wald confidence intervals
std_errors_vocation <- summary(model)$standard.errors["vocation", c("math", "science")]
ci_lower <- exp(coefficients_vocation - 1.96 * std_errors_vocation)
ci_upper <- exp(coefficients_vocation + 1.96 * std_errors_vocation)
# Combine the results into a data frame
ci_df <- data.frame(Odds_Ratio = odds_ratios, CI_Lower = ci_lower, CI_Upper = ci_upper)
ci_df
# Identify the Subject with Unexpected Coefficients
exp(summary(model)$coefficients["vocation", c("read", "write", "math", "science", "socst")])

#Q4(c)
# Create the composite score
hsb$composite <- hsb$read + hsb$write + hsb$math + hsb$science + hsb$socst
# Fit the multinomial model with the composite score
model_composite <- multinom(prog ~ gender + race + ses + schtyp + composite, data = hsb)
summary(model_composite)
# Apply stepwise AIC 
reduced_model <- step(model_composite, trace = 0)
summary(reduced_model)

#Q4(d)
# Find the most common or mean values for the predictors
most_common_gender <- names(sort(table(hsb$gender), decreasing = TRUE))[1]
most_common_race <- names(sort(table(hsb$race), decreasing = TRUE))[1]
most_common_ses <- names(sort(table(hsb$ses), decreasing = TRUE))[1]
most_common_school_type <- names(sort(table(hsb$schtyp), decreasing = TRUE))[1]
mean_composite <- mean(hsb$composite)
# Generate a sequence of composite scores over the observed range
composite_seq <- seq(min(hsb$composite), max(hsb$composite), length.out = 100)

# Create a data frame with predictors set to the most common/mean values
new_data <- data.frame(
  gender = most_common_gender,
  race = most_common_race,
  ses = most_common_ses,
  schtyp = most_common_school_type,
  composite = composite_seq
)

# Predict the probabilities for each program
predicted_probs <- predict(reduced_model, new_data, type = "probs")
# Convert predicted probabilities to a data frame for plotting
predicted_probs_df <- data.frame(composite_seq, predicted_probs)


library(ggplot2)
# Create the plot
ggplot(predicted_probs_df, aes(x = composite_seq)) +
  geom_line(aes(y = general, color = "General")) +
  geom_line(aes(y = vocation, color = "Vocational")) +
  geom_line(aes(y = academic, color = "Academic")) +
  labs(x = "Composite Score", y = "Predicted Probability", 
       title = "Predicted Probabilities of Program Choice by Composite Score") +
  scale_color_manual(values = c("General" = "blue", "Vocational" = "red", "Academic" = "green")) +
  theme_minimal()





#5(a)
library(faraway)
data(UCBAdmissions)
UCBAdmissions
# Compute odds ratios for each department
ucb_df <- data.frame(UCBAdmissions)
# Summing over the third dimension (department) to get marginal totals
ucb_marginal <- apply(UCBAdmissions, c(1, 2), sum)
ucb_marginal
OddsRatio=(1198*1278)/(1493*557) 
OddsRatio

#5(b)
data <- data.frame(UCBAdmissions)
#Admit ⊥ Gender ⊥ Dept
modelX.Y.Z. <- glm(Freq ~ Admit + Gender + Dept, family = poisson, data = data)
#XY.Z
modelXY.Z <- glm(Freq ~ Admit * Gender + Dept, family = poisson, data = data)
#Y.XZ
modelY.XZ <- glm(Freq ~ Admit * Dept + Gender, family = poisson, data = data)
#X.YZ
modelX.YZ <- glm(Freq ~ Gender * Dept + Admit, family = poisson, data = data)
#Full interaction model
modelXYZ <- glm(Freq ~ Admit * Gender * Dept, family = poisson, data = data)
#XY.YZ
modelXY.YZ <- glm(Freq ~ Admit * Gender + Gender + Dept, family = poisson, data = data)
#XZ.YZ
modelXZ.YZ <- glm(Freq ~ Admit * Dept +  Gender * Dept, family = poisson, data = data)
#XY.XZ
modelXY.XZ <- glm(Freq ~ Admit * Gender + Admit * Dept, family = poisson, data = data)
#XY.YZ.XZ
modelXY.YZ.XZ <- glm(Freq ~ Admit * Gender + Gender * Dept + Admit * Dept, family = poisson, data = data)

extract_model_stats <- function(model) {
  deviance <- deviance(model)
  pearson_chi2 <- sum(residuals(model, type = "pearson")^2)
  aic <- AIC(model)
  bic <- BIC(model)
  df_resid <- df.residual(model)
  
  return(c(Deviance = deviance, Pearson_Chi2 = pearson_chi2, AIC = aic, BIC = bic, Residual_DF = df_resid))
}

model_stats <- rbind(
  ModelX.Y.Z. = extract_model_stats(modelX.Y.Z.),
  ModelXY.Z = extract_model_stats(modelXY.Z),
  ModelY.XZ = extract_model_stats(modelY.XZ),
  ModelX.YZ = extract_model_stats(modelX.YZ),
  ModelXY.YZ = extract_model_stats(modelXY.YZ),
  ModelXZ.YZ = extract_model_stats(modelXZ.YZ),
  ModelXY.XZ = extract_model_stats(modelXY.XZ),
  ModelXY.YZ.XZ = extract_model_stats(modelXY.YZ.XZ),
  ModelXYZ = extract_model_stats(modelXYZ)
)
model_stats

#5(d)
# Extract coefficients for XZ
coef <- modelXY.YZ.XZ$coefficients["AdmitRejected:DeptC"]
# Calculate the odds ratio
OR <- exp(coef)
OR

# Extract the variance-covariance matrix
vcov_matrix <- vcov(modelXY.YZ.XZ)
se_coef <- sqrt(vcov_matrix["AdmitRejected:DeptC", "AdmitRejected:DeptC"])
# Calculate the 95% confidence interval for the log odds ratio
ci_lower_logOR <- coef - 1.96 * se_coef
ci_upper_logOR <- coef + 1.96 * se_coef

# Convert the confidence interval from log odds to odds ratio
ci_lower_OR <- exp(ci_lower_logOR)
ci_upper_OR <- exp(ci_upper_logOR)

# Display the confidence interval
paste("CI: ", ci_lower_OR, ci_upper_OR)