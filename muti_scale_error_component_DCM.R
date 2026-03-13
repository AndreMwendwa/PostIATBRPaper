# ############################################################### #
# ###################### Apollo MMNL Script ##################### #
# ############################################################### #

# Clear workspace
rm(list = ls())

# Load package
library(apollo)

# Set working directory
# setwd("./")

setwd("./")
# Initialise Apollo
apollo_initialise()

# ############################################################### #
# ######################## 1. Controls ########################## #
# ############################################################### #

apollo_control = list(
  modelName       = "MaaS_round_scale_MMNL",
  modelDescr      = "Mixed logit with round-specific scale parameters for MaaS alternatives",
  indivID         = "PersonID",
  mixing          = TRUE,
  nCores          = 63,
  outputDirectory = "output"
)

# ############################################################### #
# ######################## 2. Data ############################## #
# ############################################################### #

database = read.csv("Multimodal Options in MaaS_R1.csv", header = TRUE)

# ############################################################### #
# ################### 3. Parameters ############################# #
# ############################################################### #

apollo_beta = c(
  # Cost / time coefficients
  b_cost           = -0.15,
  b_walk           = -0.02,
  b_wait           = -0.03,
  
  # Main mode time coefficients
  b_main_bic       = -0.03,
  b_main_bs        = -0.03,
  b_main_pt        = -0.04,
  b_main_drive     = -0.01,
  b_main_ridehail  = -0.04,
  
  # Error component
  sigma_ec_BS      =  0,
  
  # Base ASC
  ASC_Walk         =  0,
  
  # Round-specific scale parameters
  mu_c1            =  1,
  mu_c2            =  1,
  mu_c3            =  1,
  
  # Random ASCs
  mu_ASC_PT        =  0,   sd_ASC_PT    =  0.5,
  mu_ASC_Drv       =  0,   sd_ASC_Drv   =  0.5,
  mu_ASC_Bic       =  0,   sd_ASC_Bic   =  0.5,
  mu_ASC_RHS       =  0,   sd_ASC_RHS   =  0.5,
  mu_ASC_BS        =  0,   sd_ASC_BS    =  0.5
)

# Fixed parameters
apollo_fixed = c("ASC_Walk", "mu_c1")

# ############################################################### #
# ###################### 4. Draws ############################### #
# ############################################################### #

apollo_draws = list(
  interDrawsType = "mlhs",
  interNDraws    = 50,
  interUnifDraws = c(),
  interNormDraws = c(
    "draws_b_main_bic",
    "draws_b_main_bs",
    "draws_b_main_pt",
    "draws_b_main_drive",
    "draws_b_main_ridehail",
    "draws_BS"
  ),
  intraDrawsType = "halton",
  intraNDraws    = 0,
  intraUnifDraws = c(),
  intraNormDraws = c()
)

# ############################################################### #
# ################### 5. Random coefficients #################### #
# ############################################################### #

apollo_randCoeff = function(apollo_beta, apollo_inputs){
  randcoeff = list()
  
  randcoeff[["ASC_PT"]]  = mu_ASC_PT  + sd_ASC_PT  * draws_b_main_pt
  randcoeff[["ASC_Drv"]] = mu_ASC_Drv + sd_ASC_Drv * draws_b_main_drive
  randcoeff[["ASC_Bic"]] = mu_ASC_Bic + sd_ASC_Bic * draws_b_main_bic
  randcoeff[["ASC_RHS"]] = mu_ASC_RHS + sd_ASC_RHS * draws_b_main_ridehail
  randcoeff[["ASC_BS"]]  = mu_ASC_BS  + sd_ASC_BS  * draws_b_main_bs
  
  randcoeff[["ec_BS"]]   = sigma_ec_BS * draws_BS
  
  return(randcoeff)
}

# ############################################################### #
# #################### 6. Validate inputs ####################### #
# ############################################################### #

apollo_inputs = apollo_validateInputs()

# ############################################################### #
# ################### 7. Probabilities ########################## #
# ############################################################### #

apollo_probabilities <- function(apollo_beta, apollo_inputs, functionality = "estimate") {
  
  apollo_attach(apollo_beta, apollo_inputs)
  on.exit(apollo_detach(apollo_beta, apollo_inputs))
  
  P <- list()
  V <- list()
  
  # ---------------- Single-mode alternatives ---------------- #
  
  V[["Drv"]] =
    ASC_Drv +
    b_main_drive * drv_time +
    b_cost * (drv_cost + park_cost)
  
  V[["PT"]] =
    ASC_PT +
    b_wait * pt_wait +
    b_main_pt * pt_time +
    b_cost * pt_cost
  
  V[["RHS"]] =
    ASC_RHS +
    b_wait * rhs_wait +
    b_main_ridehail * rhs_time +
    b_cost * rhs_cost
  
  V[["BS"]] =
    ASC_BS +
    b_main_bs * bs_time +
    b_walk * bs_walk +
    b_cost * bs_cost
  
  V[["Walk"]] =
    ASC_Walk +
    b_walk * walk_time +
    ec_BS
  
  V[["Bic"]] =
    ASC_Bic +
    b_main_bic * bic_time
  
  # ---------------- Two-mode alternatives ---------------- #
  
  V[["Bic_PT"]] =
    ASC_Bic + ASC_PT +
    b_main_bic * bicpt_bic_time +
    b_main_pt * bicpt_time +
    b_wait * bicpt_wait +
    b_cost * bicpt_cost
  
  V[["Drv_PT"]] =
    ASC_Drv + ASC_PT +
    b_main_drive * drpt_drv_time +
    b_wait * drpt_wait +
    b_main_pt * drpt_time +
    b_cost * (drpt_drv_cost + drpt_park_cost + drpt_cost)
  
  V[["PT_RHS"]] =
    ASC_PT + ASC_RHS +
    b_wait * (ptrs_wait + ptrs_rhs_wait) +
    b_main_pt * ptrs_time +
    b_main_ridehail * ptrs_rhs_time +
    b_cost * (ptrs_cost + ptrs_rhs_cost)
  
  V[["BS_PT"]] =
    ASC_BS + ASC_PT +
    b_walk * bspt_bs_walk +
    b_wait * bspt_wait +
    b_main_bs * bspt_bs_time +
    b_main_pt * bspt_time +
    b_cost * (bspt_cost + bspt_bs_cost)
  
  V[["PT_BS"]] =
    ASC_PT + ASC_BS +
    b_wait * ptbs_wait +
    b_walk * ptbs_bs_walk +
    b_main_pt * ptbs_time +
    b_main_bs * ptbs_bs_time +
    b_cost * (ptbs_bs_cost + ptbs_cost)
  
  V[["RHS_PT"]] =
    ASC_RHS + ASC_PT +
    b_main_ridehail * rhspt_rhs_time +
    b_wait * (rhspt_rhs_wait + rhspt_wait) +
    b_main_pt * rhspt_time +
    b_cost * (rhspt_cost + rhspt_rhs_cost)
  
  # ---------------- Three-mode alternatives ---------------- #
  
  V[["Bic_PT_RHS"]] =
    ASC_Bic + ASC_PT + ASC_RHS +
    b_wait * (bpr_wait + bpr_rhs_wait) +
    b_main_bic * bpr_bic_time +
    b_main_pt * bpr_time +
    b_main_ridehail * bpr_rhs_time +
    b_cost * (bpr_cost + bpr_rhs_cost)
  
  V[["RHS_PT_RHS"]] =
    ASC_RHS + ASC_PT +
    b_wait * (rpr_wait + rpr_rhs_wait + rpr_rhs2_wait) +
    b_main_ridehail * rpr_rhs_time +
    b_main_pt * rpr_time +
    b_main_ridehail * rpr_rhs2_time +
    b_cost * (rpr_cost + rpr_rhs_cost + rpr_rhs2_cost)
  
  V[["BS_PT_BS"]] =
    ASC_BS + ASC_PT +
    b_wait * bspbs_wait +
    b_walk * (bspbs_bs1_walk + bspbs_bs2_walk) +
    b_main_bs * bspbs_bs1_time +
    b_main_pt * bspbs_time +
    b_main_bs * bspbs_bs2_time +
    b_cost * (bspbs_cost + bspbs_bs1_cost + bspbs_bs2_cost)
  
  V[["BS_PT_RHS"]] =
    ASC_BS + ASC_PT + ASC_RHS +
    b_walk * bspr_bs_walk +
    b_wait * (bspr_wait + bspr_rhs_wait) +
    b_main_bs * bspr_bs_time +
    b_main_pt * bspr_time +
    b_main_ridehail * bspr_rhs_time +
    b_cost * (bspr_cost + bspr_rhs_cost + bspr_bs_cost)
  
  V[["Drv_PT_BS"]] =
    ASC_Drv + ASC_PT + ASC_BS +
    b_main_drive * dpbs_drv_time +
    b_walk * dpbs_bs_walk +
    b_wait * dpbs_wait +
    b_main_pt * dpbs_time +
    b_main_bs * dpbs_bs_time +
    b_cost * (dpbs_bs_cost + dpbs_drv_cost + dpbs_drv_park + dpbs_cost)
  
  V[["RHS_PT_BS"]] =
    ASC_RHS + ASC_PT + ASC_BS +
    b_walk * rpbs_bs_walk +
    b_wait * (rpbs_wait + rpbs_rhs_wait) +
    b_main_pt * rpbs_time +
    b_main_ridehail * rpbs_rhs_time +
    b_main_bs * rpbs_bs_time +
    b_cost * (rpbs_bs_cost + rpbs_cost + rpbs_rhs_cost)
  
  V[["Bic_PT_BS"]] =
    ASC_Bic + ASC_PT + ASC_BS +
    b_walk * bpbs_bs_walk +
    b_wait * bpbs_wait +
    b_main_pt * bpbs_time +
    b_main_bic * bpbs_bic_time +
    b_main_bs * bpbs_bs_time +
    b_cost * (bpbs_cost + bpbs_bs_cost)
  
  V[["Drv_PT_RHS"]] =
    ASC_Drv + ASC_PT + ASC_RHS +
    b_main_drive * dpr_drv_time +
    b_wait * (dpr_wait + dpr_rhs_wait) +
    b_main_pt * dpr_time +
    b_main_ridehail * dpr_rhs_time +
    b_cost * (dpr_cost + dpr_rhs_cost + dpr_drv_cost + dpr_park_cost)
  
  # ---------------- Shared alternative / availability definitions ---------------- #
  
  alternatives = c(
    Bic_PT_RHS = 0,
    Bic_PT = 1,
    PT = 2,
    Drv_PT = 3,
    Drv_PT_RHS = 4,
    PT_RHS = 5,
    Drv = 6,
    RHS_PT = 7,
    Walk = 8,
    RHS_PT_RHS = 9,
    BS_PT_BS = 10,
    BS_PT_RHS = 11,
    PT_BS = 12,
    Bic = 13,
    RHS = 14,
    Drv_PT_BS = 15,
    RHS_PT_BS = 16,
    BS_PT = 17,
    Bic_PT_BS = 18,
    BS = 19
  )
  
  avail_list = list(
    Bic_PT_RHS = av_Bic_PT_RHS,
    Bic_PT = av_Bic_PT,
    PT = av_PT,
    Drv_PT = av_Drv_PT,
    Drv_PT_RHS = av_Drv_PT_RHS,
    PT_RHS = av_PT_RHS,
    Drv = av_Drv,
    RHS_PT = av_RHS_PT,
    Walk = av_Walk,
    RHS_PT_RHS = av_RHS_PT_RHS,
    BS_PT_BS = av_BS_PT_BS,
    BS_PT_RHS = av_BS_PT_RHS,
    PT_BS = av_PT_BS,
    Bic = av_Bic,
    RHS = av_RHS,
    Drv_PT_BS = av_Drv_PT_BS,
    RHS_PT_BS = av_RHS_PT_BS,
    BS_PT = av_BS_PT,
    Bic_PT_BS = av_Bic_PT_BS,
    BS = av_BS
  )
  
  make_utilities <- function(mu_scale, V){
    list(
      Bic_PT_RHS = mu_scale * V[["Bic_PT_RHS"]],
      Bic_PT = mu_scale * V[["Bic_PT"]],
      PT = mu_scale * V[["PT"]],
      Drv_PT = mu_scale * V[["Drv_PT"]],
      Drv_PT_RHS = mu_scale * V[["Drv_PT_RHS"]],
      PT_RHS = mu_scale * V[["PT_RHS"]],
      Drv = mu_scale * V[["Drv"]],
      RHS_PT = mu_scale * V[["RHS_PT"]],
      Walk = mu_scale * V[["Walk"]],
      RHS_PT_RHS = mu_scale * V[["RHS_PT_RHS"]],
      BS_PT_BS = mu_scale * V[["BS_PT_BS"]],
      BS_PT_RHS = mu_scale * V[["BS_PT_RHS"]],
      PT_BS = mu_scale * V[["PT_BS"]],
      Bic = mu_scale * V[["Bic"]],
      RHS = mu_scale * V[["RHS"]],
      Drv_PT_BS = mu_scale * V[["Drv_PT_BS"]],
      RHS_PT_BS = mu_scale * V[["RHS_PT_BS"]],
      BS_PT = mu_scale * V[["BS_PT"]],
      Bic_PT_BS = mu_scale * V[["Bic_PT_BS"]],
      BS = mu_scale * V[["BS"]]
    )
  }
  
  # ---------------- Round 1 ---------------- #
  
  mnl_settings_c1 <- list(
    alternatives = alternatives,
    avail        = avail_list,
    choiceVar    = choice,
    utilities    = make_utilities(mu_c1, V),
    rows         = (rounds == 1)
  )
  
  # ---------------- Round 2 ---------------- #
  
  mnl_settings_c2 <- list(
    alternatives = alternatives,
    avail        = avail_list,
    choiceVar    = choice,
    utilities    = make_utilities(mu_c2, V),
    rows         = (rounds == 2)
  )
  
  # ---------------- Round 3 ---------------- #
  
  mnl_settings_c3 <- list(
    alternatives = alternatives,
    avail        = avail_list,
    choiceVar    = choice,
    utilities    = make_utilities(mu_c3, V),
    rows         = (rounds == 3)
  )
  
  # ---------------- MNL probabilities ---------------- #
  
  P[["c1"]] = apollo_mnl(mnl_settings_c1, functionality)
  P[["c2"]] = apollo_mnl(mnl_settings_c2, functionality)
  P[["c3"]] = apollo_mnl(mnl_settings_c3, functionality)
  
  # Combine across rounds
  P = apollo_combineModels(P, apollo_inputs, functionality)
  
  # Panel product across observations for each individual
  P = apollo_panelProd(P, apollo_inputs, functionality)
  
  # Average across inter-individual draws
  P = apollo_avgInterDraws(P, apollo_inputs, functionality)
  
  # Final preparation
  P = apollo_prepareProb(P, apollo_inputs, functionality)
  
  return(P)
}

# ############################################################### #
# ###################### 8. Estimation ########################## #
# ############################################################### #

model = apollo_estimate(
  apollo_beta,
  apollo_fixed,
  apollo_probabilities,
  apollo_inputs
)

# ############################################################### #
# ######################## 9. Output ############################ #
# ############################################################### #

apollo_modelOutput(model)
apollo_saveOutput(model)








