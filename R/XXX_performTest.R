

perform_test.default <- function(x, y, alternative = c("two.sided", "greater", "less"), delta = ifelse(var.test, 1, 0),
                         test.name, method = c("asymptotic", "permutation", "randomization"), n.rep = ifelse(method == "randomization", 10000, NULL),
                         na.rm = FALSE, var.test = FALSE) {

  if (!(test.name %in% c("HL11", "HL12", "HL21", "HL22", "MED1", "MED2"))) {
    stop(paste("Test", test.name, "not recognized. Please use one of the following: 'HL11', 'HL12', 'HL21', 'HL22', 'MED1', 'MED2'."))
  }

  alternative <- match.arg(alternative)
  method <- match.arg(method)

  if (test.name == "HL11") {
    hl1_test(x = x, y = y, alternative = alternative, delta = delta, method = method, scale = "S1", n.rep = n.rep, na.rm = na.rm, var.test = var.test)
  } else if (test.name == "HL12") {
    hl1_test(x = x, y = y, alternative = alternative, delta = delta, method = method, scale = "S2", n.rep = n.rep, na.rm = na.rm, var.test = var.test)
  } else if (test.name == "HL21") {
    hl2_test(x = x, y = y, alternative = alternative, delta = delta, method = method, scale = "S1", n.rep = n.rep, na.rm = na.rm, var.test = var.test)
  } else if (test.name == "HL22") {
    hl2_test(x = x, y = y, alternative = alternative, delta = delta, method = method, scale = "S2", n.rep = n.rep, na.rm = na.rm, var.test = var.test)
  } else if (test.name == "MED1") {
    med_test(x = x, y = y, alternative = alternative, delta = delta, method = method, scale = "S3", n.rep = n.rep, na.rm = na.rm, var.test = var.test)
  } else if (test.name == "MED2") {
    med_test(x = x, y = y, alternative = alternative, delta = delta, method = method, scale = "S4", n.rep = n.rep, na.rm = na.rm, var.test = var.test)
  }


}
