#' For quick checking the version of packages required
library(devtools)

# mlr3 and others
package_info('mlr3', dependencies = F) # 0.15.0
package_info('mlr3pipelines', dependencies = F) # 0.4.3
package_info('mlr3tuning', dependencies = F) # 0.18.0
package_info('paradox', dependencies = F) # 0.11.1
package_info('mlr3mbo', dependencies = F) # 0.2.0@e1deed
package_info('DiceKriging', dependencies = F) # 1.6.0
package_info('mlr3filters', dependencies = F) # 0.7.1@c684ecf
#' fixing `find_correlation` filter that couldn't be used with survival tasks
package_info('mlr3proba', dependencies = F) # 0.5.2@e263908
#' lots of fixes regarding measures (RCLL, ERV measures, etc.)
package_info('distr6', dependencies = F) # 1.6.15
package_info('stabm', dependencies = F) # 1.2.2
package_info('Matrix', dependencies = F) # 1.5-4
package_info('tictoc', dependencies = F) # 1.1

# model-related
package_info('mlr3extralearners', dependencies = F) # 0.6.1
package_info('glmnet', dependencies = F) # 4.1-7
package_info('ranger', dependencies = F) # 0.14.1
package_info('CoxBoost', dependencies = F) # binderh/CoxBoost@1dc47d7
package_info('aorsf', dependencies = F) # 0.0.7
package_info('xgboost', dependencies = F) # 1.7.5.1
