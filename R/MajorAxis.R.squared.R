# ##R-squared calculations for major axis analysis##
#
# ##Calculate the multivariate sum of squared differences of PC scores from the mean PC scores
# ##Calculate the multivariate sum of squared differences of PC scores from the fitted value for that specimens' PC1 score
# ##Calculate r.squared by dividing top from bottom and subtracting this from 1
#
# X <- Major_Axis_lm
#
#
# ##Loading in temp values
# temp_PC1_values <- Combined_Morphospace_Clades3$PCvalues$`Terrestrial Notosuchia`[,"PC1"]
# temp_PC2_values <- Combined_Morphospace_Clades3$PCvalues$`Terrestrial Notosuchia`[,"PC2"]
# temp_slope <- Clades3.MajorAxis$Model_list$`Terrestrial Notosuchia`[[1]]$coefficients[2]
# temp_intercept <- Clades3.MajorAxis$Model_list$`Terrestrial Notosuchia`[[1]]$coefficients[1]
# temp_fit_values <- (temp_PC1_values * temp_slope) + temp_intercept
#
# ##Calculate the multivariate sum of squared differences of PC scores from the mean PC scores
# temp_SS_tot <- sum((temp_PC2_values - mean(temp_PC2_values))^2)
#
# ##Calculate the multivariate sum of squared differences of PC scores from the fitted value for that specimens' PC1 score
# temp_SS_res <- sum((temp_PC2_values - temp_fit_values)^2)
#
# ##Calculate r.squared by dividing top from bottom and subtracting this from 1
# temp_r.squared <- 1 - (temp_SS_res/temp_SS_tot)
#
# ###
#
# Clades3.MajorAxis$Model_list$`Terrestrial Notosuchia`[[1]]$coefficients[1]
# Clades3.MajorAxis$Model_list$`Terrestrial Notosuchia`[[1]]$coefficients[2]
#
# summary(lm(Combined_Morphospace_Clades3$PCvalues$`Terrestrial Notosuchia`[,"PC2"] ~ Combined_Morphospace_Clades3$PCvalues$`Terrestrial Notosuchia`[,"PC1"]))
#
