
##Preparatory steps##
setwd()

install.packages(c("tidyverse", "Amelia", "stats", "dplyr", "survival", "survminer", "broom"))

library(tidyverse)
library(Amelia)
library(stats)
library(dplyr)
library(survival)
library(survminer)
library(broom)

# Read the dataset (ensure the file is in the working directory or provide the full path)
bmi <- read_dta("bmi.dta")

# View the structure of the dataset
str(bmi)

# View the first few rows of the dataset
head(bmi)


##Checking the model assumptions

#3. re-fit the model and obtain the predicted values
ls()
str(bmi)
bmi_data <- bmi


# Fit the linear model
model <- lm(bmi ~ age, data = bmi_data)

# Add predictions to the dataset
bmi_data$yhat <- predict(model)

#4. Calculate the studentised residuals

bmi_data$resid <- rstudent(model)
hist(bmi_data$resid, main = "Histogram of Residuals", xlab = "Residuals", col = "blue", border = "white")
 
#5. 5.	Draw a quantile-quantile (QQ) plot of the residuals, using the command:

qqnorm(bmi_data$resid)
qqline(bmi_data$resid, col = "red")

#6. 6.	Graph the residuals against the fitted values
ggplot(bmi_data, aes(x = yhat, y = resid)) + geom_point() +
  labs(title = "Residuals vs Fitted Values", x = "Fitted values", y = "Residuals") +
  theme_minimal()

#8.8.	We can perform a hypothesis test for heteroskedasticity
library(lmtest)
bptest(model)

##Influential data points##
#9. high leverage or influence

# Leverage plot
library(ggplot2)
library(car)
bmi_data$leverage <- hatvalues(model) # Calculate leverage (hat) values 
bmi_data$residuals_raw <- residuals(model) # Extract raw residuals
rss <- sum(bmi_data$residuals_raw^2) # Calculate residual sum of squares

bmi_data$norm_resid <- (bmi_data$residuals_raw) / sqrt(rss) # Calculate normalized residuals squared
bmi_data$norm_resid_sq <- bmi_data$norm_resid^2 # Calculate normalized residuals squared
ggplot(bmi_data, aes(x = norm_resid_sq, y = leverage)) + geom_point() +
  labs(title = "Leverage vs. Squared Standardized Residuals", x = "Normalized Residuals Squared",
       y = "Leverage") + theme_minimal()
# R Alternative (influence plot showing leverage and Cook's distance) 
influencePlot(model, main = "Influence Plot", id.method = "identify", id.n = 5)

#10. calculate ??????? (dfbetas), to assess the influence of individual observations on the regression coefficient for age
dfbetas_values <- dfbetas(model)


#12.	Finally, use the following commands to graph the data
# Calculate absolute DFBETA values for 'age' predictor
bmi_data$abs_dfbeta <- abs(bmi_data$dfbeta_1)

# Create the plot
ggplot(bmi_data, aes(x = age, y = bmi)) + 
  geom_point(aes(size = abs_dfbeta), shape = 1) +  # Circles with sizes proportional to abs_dfbeta
  geom_line(aes(y = yhat), color = "red") +       # Line for fitted values (yhat)
  labs(title = "Relationship Between BMI & Age with DFBETA Weighting",
       x = "Age", 
       y = "BMI") + 
  theme_minimal() + 
  scale_size_continuous(range = c(1, 10))  # Change size of circles


###PART 2 - ROBUST FIT
# Preparatory steps

install.packages("dplyr")   # Install dplyr for data manipulation
install.packages("sandwich")  # Install sandwich for robust standard errors
install.packages("haven")  # If not already installed


library(haven)  # Load the haven package


#
library(foreign)
library(dplyr)
library(ggplot2)
library(sandwich)
library(lmtest)
#
# Load the dataset

vitd_data <- read_dta("vitd.dta")  # Adjust the path to the location of vitd.dta



###Data exploration##
# Summarize the vitamin D measurements
summary(vitd_data$vitd)

# Plot histogram
library(ggplot2)
ggplot(vitd_data, aes(x = vitd)) +
  geom_histogram(binwidth = 5, fill = "blue", color = "black") +
  labs(title = "Histogram of Vitamin D Measurements", x = "Vitamin D Level", y = "Frequency") +
  theme_minimal()

#3. Variable t, representing the time of the year the vitamin D measurement was taken 
# Install dplyr if not installed
install.packages("dplyr")

# Load dplyr package
library(dplyr)

vitd_data <- vitd_data %>%
  mutate(year = as.integer(year), # Ensure year is integer
         date_test = as.Date(date_test, format = "%Y-%m-%d"), # Parse date if needed
         t = as.numeric((date_test - as.Date(paste0("01jan", year), format = "%d%b%Y")) / 365.25))

#4. 4.	Draw a scatter plot of vitamin D and time.
install.packages("ggplot2")
library(ggplot2)  # Load ggplot2 package

ggplot(vitd_data, aes(x = t, y = vitd)) +
  geom_point(color = "blue") +
  labs(title = "Scatter Plot of Vitamin D vs. Time", x = "Time (t)", y = "Vitamin D") +
  theme_minimal()


#5.	Fit quadratic regression model of vitamin D and time
# Install and load necessary packages
install.packages("dplyr")
install.packages("ggplot2")
library(dplyr)
library(ggplot2)
library(haven)


# Fit quadratic regression model
fit <- lm(vitd ~ t + I(t^2), data = vitd_data)

# Add predicted values to the dataset
vitd_data$vitd_pred <- predict(fit)

# Plot scatter and fitted quadratic regression line
ggplot(vitd_data, aes(x = t)) +
  geom_point(aes(y = vitd), color = "blue") +
  geom_line(aes(y = vitd_pred), color = "red") +
  labs(title = "Scatter of Vitamin D and Fitted Values", x = "Time (t)", y = "Vitamin D") +
  theme_minimal()
summary(fit)


##More robust inference##
#6	Create a residual-vs-fitted value plot for the previous regression model. Does the assumption of constant variance appear to hold?
# Residuals vs. Fitted values
ggplot(data.frame(fitted = fitted(fit), residuals = residuals(fit)), aes(x = fitted, y = residuals)) +
  geom_point(color = "blue") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "Residuals vs. Fitted Values", x = "Fitted Values", y = "Residuals") +
  theme_minimal()

#7.	Refit the model, calculating robust standard errors

# Install lmtest package if not already installed
if (!requireNamespace("lmtest", quietly = TRUE)) {
  install.packages("lmtest")
}

# Load the package
library(lmtest)


# Compute robust standard errors
robust_se <- coeftest(fit, vcov = vcovHC(fit, type = "HC1"))

# Compare regular and robust standard errors
print("Regular Standard Errors:")
print(summary(fit)$coefficients)

print("Robust Standard Errors:")
print(robust_se)

