createMs = function(preds) {
  
  # create all combinations models
  if (length(preds) == 4) {
    variables = c(preds[1], preds[2], preds[3], preds[4], 
                  sprintf("%s:%s", preds[1], preds[2]), 
                  sprintf("%s:%s", preds[1], preds[3]), 
                  sprintf("%s:%s", preds[1], preds[4]),
                  sprintf("%s:%s", preds[2], preds[3]), 
                  sprintf("%s:%s", preds[2], preds[4]), 
                  sprintf("%s:%s", preds[3], preds[4]), 
                  sprintf("%s:%s:%s", preds[1], preds[2], preds[3]),
                  sprintf("%s:%s:%s", preds[1], preds[2], preds[4]), 
                  sprintf("%s:%s:%s", preds[2], preds[3], preds[4]), 
                  sprintf("%s:%s:%s", preds[1], preds[3], preds[4]))
  } else if (length(preds) == 3) {
    variables = c(preds[1], preds[2], preds[3], 
                  sprintf("%s:%s", preds[1], preds[2]), 
                  sprintf("%s:%s", preds[1], preds[3]), 
                  sprintf("%s:%s", preds[2], preds[3]), 
                  sprintf("%s:%s:%s", preds[1], preds[2], preds[3]))
  }
  
  formulas = c()
  for (i in seq_along(variables)) {
    tmp = combn(variables, i)
    tmp = apply(tmp, 2, paste, collapse=" + ")
    formulas = c(formulas, tmp)
  }
  
  # get rid of impossible models
  fixed = c()
  for (f in formulas) {
    # start with assumption that it is okay
    use = T
    # extract predictors
    preds = str_split_1(f, pattern = " \\+ ")
    for (p in preds) {
      # check if it's an interaction term
      if (grepl(":", p)) {
        # split into the predictors included in the interaction
        intpreds = str_split_1(p, ":")
        # check for each of the predictors, if the main effect is included
        for (i in intpreds) {
          if (!any(preds == i)) {
            use = F
            next
          }
        }
        # check for three-way interactions if two-ways are included
        if (length(intpreds) > 2) {
          if (!any(preds == sprintf('%s:%s', intpreds[1], intpreds[2])) |
              !any(preds == sprintf('%s:%s', intpreds[1], intpreds[3])) |
              !any(preds == sprintf('%s:%s', intpreds[2], intpreds[3]))) {
            use = F
            next
          }
        }
      }
    }
    if (use) fixed = c(fixed, f)
  }
  
  return(fixed)
  
}
