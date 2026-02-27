

# Set the working directory
library(apollo)
setwd("./")

####################################################
####################################################
################### MNL Example ####################
####################################################
####################################################

################### Step 1 #########################
# Clear the work space
rm(list = ls())
# cdes
# Initialize the MNL Code 
apollo_initialise()

################### Step 2 #########################
# Set the core controls   
### Set core controls
apollo_control = list(
  modelDescr      = "Round 1 only",
  indivID         = "PersonID",  
  nCores          = 63,
  outputDirectory = "output"
)

# Read in the mode choice data
database= read.csv("Multimodal Options in MaaS_R1.csv", header = TRUE)

# Restricting to round 1
database_1 <- database[database$rounds == 1, ]
database <- database_1

################### Step 3 #########################
# Define the model parameters and set the starting value
### -- 1) ???????????? -- ###
apollo_beta = c(
  # -- ??????????????? (log-normal) -- #
  b_cost       = -0.15,    # ????????????????????????
  
  # -- ???????????? (normal) -- #
  b_walk       = -0.02,    # ???????????????
  b_wait       = -0.03,    # ???????????????
  
  # -- ???????????? (normal) -- #
  b_main_bic       = -0.03,
  b_main_bs        = -0.03,
  b_main_pt        = -0.04,
  b_main_drive     = -0.01,
  b_main_ridehail  = -0.04,
  b_tr             =  0,
  
  # -- ?????????????????? -- #
  sigma_ec_BS      =  0,
  
  # -- ???????????? (??? Walk ?????? 0) -- #
  ASC_Walk         =  0,
  
  # -- ?????????????????????/????????? -- #
  mu_ASC_PT        =  0,   sd_ASC_PT    =  0.5,
  mu_ASC_Drv       =  0,   sd_ASC_Drv   =  0.5,
  mu_ASC_Bic       =  0,   sd_ASC_Bic   =  0.5,
  mu_ASC_RHS       =  0,   sd_ASC_RHS   =  0.5,
  mu_ASC_BS        =  0,   sd_ASC_BS    =  0.5
)


# Define the model parameters that will be fixed to the starting values
apollo_fixed = c("ASC_Walk")


apollo_draws = list(
  interDrawsType = "mlhs",
  interNDraws    = 50,
  interUnifDraws = c(),
  interNormDraws  = c(
    
    
    
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


### -- 3) ?????????????????? -- ###
apollo_randCoeff = function(apollo_beta, apollo_inputs){
  randcoeff = list()
  # log-normal ??????:??????
  
  # normal ??????
  
  randcoeff[["ASC_PT"]] = mu_ASC_PT + sd_ASC_PT * draws_b_main_pt
  randcoeff[["ASC_Drv"]] = mu_ASC_Drv + sd_ASC_Drv * draws_b_main_drive
  randcoeff[["ASC_Bic"]] = mu_ASC_Bic + sd_ASC_Bic * draws_b_main_bic
  randcoeff[["ASC_RHS"]] = mu_ASC_RHS + sd_ASC_RHS * draws_b_main_ridehail
  randcoeff[["ASC_BS"]] = mu_ASC_BS + sd_ASC_BS * draws_b_main_bs
  
  
  # ??????????????????
  
  randcoeff[["ec_BS"]]           = sigma_ec_BS  * draws_BS
  
  
  return(randcoeff)}



################### Step 4 #########################
# Validate the input data to ensure that there are no errors
apollo_inputs = apollo_validateInputs()

################### Step 5 #########################
# Define the model and the likelihood function
apollo_probabilities <- function(apollo_beta, apollo_inputs, functionality = "estimate") {
  
  ### Attach inputs and detach after function exit
  apollo_attach(apollo_beta, apollo_inputs)
  on.exit(apollo_detach(apollo_beta, apollo_inputs))
  
  ### Create list of probabilities P
  P <- list()
  
  #
  #-- ?????????????????? V ??????,?????????????????? Access/Main/main ????????????????????? --#
  V <- list()
  
  #-- 1. ????????? Drv --#
  V[['Drv']] =
    ASC_Drv +
    b_main_drive * drv_time +
    b_cost       * (drv_cost + park_cost)
  # no EC here
  
  #-- 2. ?????????/?????? PT --#
  V[['PT']] =
    ASC_PT +
    b_wait * pt_wait +
    b_main_pt   * pt_time +
    b_cost      * pt_cost 
  
  #-- 3. ???????????? RHS --#
  V[['RHS']] =
    ASC_RHS +
    b_wait       * rhs_wait +
    b_main_ridehail   * rhs_time +
    b_cost            * rhs_cost 
  
  #-- 4. ??????????????? BS --#
  V[['BS']] =
    ASC_BS +
    b_main_bs   * bs_time +
    b_walk * bs_walk +
    b_cost      * bs_cost 
  
  
  #-- 5. ????????? Walk --#
  V[['Walk']] =
    ASC_Walk +
    b_walk * walk_time+ec_BS
  # no EC
  
  #-- 6. ????????? Bic --#
  V[['Bic']] =
    ASC_Bic +
    b_main_bic * bic_time
  # no EC
  
  #-- 7. ??????????????? Bic_PT --#
  V[['Bic_PT']] =
    ASC_Bic +ASC_PT +
    b_main_bic     * bicpt_bic_time +
    b_main_pt      * bicpt_time +
    b_wait    * bicpt_wait +
    b_cost         * bicpt_cost 
    + b_tr * 1
  
  
  
  #-- 8. ??????????????? Drv_PT --#
  V[['Drv_PT']] =
    ASC_Drv +ASC_PT+
    b_main_drive     * drpt_drv_time +
    b_wait      * drpt_wait +
    b_main_pt        * drpt_time +
    b_cost           * (drpt_drv_cost + drpt_park_cost + drpt_cost) 
  + b_tr * 1
  
  #-- 9. ?????????????????? PT_RHS --#
  V[['PT_RHS']] =
    ASC_PT +ASC_RHS +
    b_wait        * (ptrs_wait + ptrs_rhs_wait) +
    b_main_pt          * ptrs_time +
    b_main_ridehail    * ptrs_rhs_time +
    b_cost             * (ptrs_cost + ptrs_rhs_cost)
  
  #-- 10. ????????????????????? BS_PT --#
  V[['BS_PT']] =
    ASC_BS +ASC_PT +
    b_walk       * bspt_bs_walk +b_wait* bspt_wait+
    b_main_bs         * bspt_bs_time +
    b_main_pt         * bspt_time +
    b_cost            * (bspt_cost + bspt_bs_cost) 
  + b_tr * 1
  
  #-- 11. ????????????????????? PT_BS --#
  V[['PT_BS']] =
    ASC_PT +ASC_BS +
    b_wait   * (ptbs_wait ) +  b_walk       * (ptbs_bs_walk)+
    
    b_main_pt         * ptbs_time +
    b_main_bs         * ptbs_bs_time +
    b_cost            * (ptbs_bs_cost + ptbs_cost) 
  + b_tr * 1
  
  #-- 12. ?????????????????? RHS_PT --#
  V[['RHS_PT']] =
    ASC_RHS +ASC_PT +
    b_main_ridehail  * rhspt_rhs_time +
    b_wait      * (rhspt_rhs_wait + rhspt_wait) +
    b_main_pt        * rhspt_time +
    b_cost           * (rhspt_cost + rhspt_rhs_cost)
  + b_tr * 1
  
  #-- 13. ??????????????????????????? Bic_PT_RHS --#
  V[['Bic_PT_RHS']] =
    ASC_Bic +ASC_PT +ASC_RHS +
    b_wait       * (bpr_wait + bpr_rhs_wait) +
    b_main_bic        * bpr_bic_time +
    b_main_pt         * bpr_time +
    b_main_ridehail   * bpr_rhs_time +
    b_cost            * (bpr_cost + bpr_rhs_cost) 
  + b_tr * 2
  
  #-- 14. ?????????????????????????????? RHS_PT_RHS --#
  V[['RHS_PT_RHS']] =
    ASC_RHS +ASC_PT +
    b_wait       * (rpr_wait + rpr_rhs_wait + rpr_rhs2_wait) +
    b_main_ridehail   * rpr_rhs_time +
    b_main_pt         * rpr_time +
    b_main_ridehail   * rpr_rhs2_time +
    b_cost            * (rpr_cost + rpr_rhs_cost + rpr_rhs2_cost)
  + b_tr * 2
  
  #-- 15. BS???PT???BS BS_PT_BS --#
  V[['BS_PT_BS']] =
    ASC_BS + ASC_PT+
    b_wait     * (bspbs_wait ) +
    b_walk     * (bspbs_bs1_walk + bspbs_bs2_walk) +
    
    b_main_bs       * bspbs_bs1_time +
    b_main_pt       * bspbs_time +
    b_main_bs       * bspbs_bs2_time +
    b_cost          * (bspbs_cost + bspbs_bs1_cost + bspbs_bs2_cost) 
  + b_tr * 2
  
  #-- 16. BS???PT???RHS BS_PT_RHS --#
  V[['BS_PT_RHS']] =
    ASC_BS +ASC_PT +ASC_RHS +
    b_walk        * (bspr_bs_walk) +
    b_wait        * (bspr_wait + bspr_rhs_wait ) +
    
    b_main_bs          * bspr_bs_time +
    b_main_pt          * bspr_time +
    b_main_ridehail    * bspr_rhs_time +
    b_cost             * (bspr_cost + bspr_rhs_cost + bspr_bs_cost) 
  + b_tr * 2
  
  #-- 17. Drv???PT???BS Drv_PT_BS --#
  V[['Drv_PT_BS']] =
    ASC_Drv +ASC_PT +ASC_BS +
    b_main_drive       * dpbs_drv_time +
    
    b_walk        * (dpbs_bs_walk) +
    b_wait        * (dpbs_wait) +
    
    b_main_pt          * dpbs_time +
    b_main_bs          * dpbs_bs_time +
    b_cost             * (dpbs_bs_cost + dpbs_drv_cost + dpbs_drv_park + dpbs_cost)
  + b_tr * 2
  
  #-- 18. RHS???PT???BS RHS_PT_BS --#
  V[['RHS_PT_BS']] =
    ASC_RHS +ASC_PT +ASC_BS +
    b_walk        * (rpbs_bs_walk) +
    b_wait        * (rpbs_wait + rpbs_rhs_wait) +
    b_main_pt          * rpbs_time +
    b_main_ridehail    * rpbs_rhs_time +
    b_main_bs          * rpbs_bs_time +
    b_cost             * (rpbs_bs_cost + rpbs_cost + rpbs_rhs_cost) 
  + b_tr * 2
  
  #-- 19. Bic???PT???BS Bic_PT_BS --#
  V[['Bic_PT_BS']] =
    ASC_Bic +ASC_PT+ASC_BS+
    b_walk        * (bpbs_bs_walk ) +
    b_wait        * (bpbs_wait) +
    b_main_pt          * bpbs_time +
    b_main_bic         * bpbs_bic_time +
    b_main_bs          * bpbs_bs_time +
    b_cost             * (bpbs_cost + bpbs_bs_cost)
  + b_tr * 2
  
  #-- 20. Drv???PT???RHS Drv_PT_RHS --#
  V[['Drv_PT_RHS']] =
    ASC_Drv +ASC_PT+ASC_RHS+
    b_main_drive       * dpr_drv_time +
    b_wait        * (dpr_wait + dpr_rhs_wait) +
    b_main_pt          * dpr_time +
    b_main_ridehail    * dpr_rhs_time +
    b_cost             * (dpr_cost + dpr_rhs_cost + dpr_drv_cost + dpr_park_cost)
  + b_tr * 2
  
  
  ### Define settings for MNL model component
  mnl_settings <- list(
    alternatives = c(Bic_PT_RHS= 0,
                     Bic_PT= 1,
                     PT= 2,
                     Drv_PT= 3,
                     Drv_PT_RHS= 4,
                     PT_RHS= 5,
                     Drv= 6,
                     RHS_PT= 7,
                     Walk= 8,
                     RHS_PT_RHS= 9,
                     BS_PT_BS= 10,
                     BS_PT_RHS= 11,
                     PT_BS= 12,
                     Bic= 13,
                     RHS= 14,
                     Drv_PT_BS= 15,
                     RHS_PT_BS= 16,
                     BS_PT= 17,
                     Bic_PT_BS= 18,
                     BS= 19), 
    
    avail = list(    Bic_PT_RHS = av_Bic_PT_RHS,
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
                     BS = av_BS), 
    choiceVar = choice,
    utilities = V
  )
  
  ### Compute probabilities using MNL model
  P[['model']] = apollo_mnl(mnl_settings, functionality)
  
  ### Take product across observation for same individual
  P = apollo_panelProd(P, apollo_inputs, functionality)
  
  ### Average across inter-individual draws
  P = apollo_avgInterDraws(P, apollo_inputs, functionality)
  
  ### Prepare and return outputs of function
  P = apollo_prepareProb(P, apollo_inputs, functionality)
  return(P)
}

# ################################################################# #
#### MODEL ESTIMATION                                            ####
# ################################################################# #

model = apollo_estimate(apollo_beta, apollo_fixed,apollo_probabilities, apollo_inputs)

# ################################################################# #
#### MODEL OUTPUTS                                               ####
# ################################################################# #

# ----------------------------------------------------------------- #
#---- FORMATTED OUTPUT (TO SCREEN)                               ----
# ----------------------------------------------------------------- #

apollo_modelOutput(model)

# ----------------------------------------------------------------- #
#---- FORMATTED OUTPUT (TO FILE, using model name)               ----
# ----------------------------------------------------------------- #

apollo_saveOutput(model)



