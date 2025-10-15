library(readr)

# ==== Load data ====
df <- read_csv("log_lifespan_predictions.csv")

# ==== Filter only Test dataset ====
df_test <- df[df$Dataset == "Test", ]

# ==== Define reusable R² function ====
calc_r2 <- function(y, yhat) {
  SSres <- sum((y - yhat)^2)
  SStot <- sum((y - mean(y))^2)
  return(1 - SSres / SStot)
}

# ==== Calculate R² for Model Averaging predictions ====
r2_value <- calc_r2(
  df_test$Known_Log_Lifespan,
  df_test$Predicted_Log_Lifespan_Model_Averaging
)

print(r2_value)
