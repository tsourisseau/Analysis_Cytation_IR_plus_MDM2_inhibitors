#Analysis of Imed's results "RT_Cytation_20250213"
#Reading 4 days after IR

#install.packages("drc")
#install.packages("ggplot2")

library(drc)
library(ggplot2)
library(dplyr)

setwd("//cifs/partage_sarcome/Equipe Italiano/Epigénétique_Tony/MDM2_IR/RT_Cytation_20250213_Imed_AnalysisTS")



#initialise experiment parameters
dose <- list("0Gy", "2Gy", "4Gy", "6Gy", "8Gy")
ttt <- list("KTX169", "BI907828")
Conc_KTX169 <- rep(x = c(0, 1.024E-3, 5.12E-3, 0.0256, 0.128, 0.64, 3.2, 16, 80, 400, 0),
                   times = 4)
Conc_BI907828 <- rep(x = c(0, 0, 0.00128, 0.0064, 0.032, 0.16, 0.8, 4, 20, 100, 500),
                     times = 4)
cell_lines <- list("IB114", "IB115", "IB128")
replicates <- 4  # 4 replicates per cell line


#read the 5 .txt files and save as variables
#The 5 .xls files (5 x 384-well plates) were saved as .txt files in one directory
List_Data_files <- list.files("\\\\cifs/partage_sarcome/Equipe Italiano/Epigénétique_Tony/MDM2_IR/RT_Cytation_20250213_Imed_AnalysisTS/RT_Cytation_20250213_Imed_AnalysisTS_RawData")
for (i in c(1:length(List_Data_files))) {
  #fetch data from file
  Data_nGy <- read.table(file = paste0("RT_Cytation_20250213_Imed_AnalysisTS_RawData/", 
                                       List_Data_files[i]),
                         sep = "\t",
                         na.strings = "NA")
  #generate data name
  Data_nGy_name <- paste0("Data_", dose[[i]])
  #assign data
  assign(x = Data_nGy_name, value = Data_nGy)
}

#compile the 5 data tables in one list of tables
RawData <- list("0Gy" = Data_0Gy,
                "2Gy" = Data_2Gy,
                "4Gy" = Data_4Gy,
                "6Gy" = Data_6Gy,
                "8Gy" = Data_8Gy)
#save RawData file to export for synergyfinder
#save(RawData, file = "RawData_RT_Cytation_20250213_Imed.RData")

#remove raw data
i_dose <- NULL
data_to_remove <- c()
for (i_dose in c(1:length(dose))) {
  data_to_remove <- c(data_to_remove, paste0("Data_", dose[[i_dose]]))
}
rm(list = data_to_remove)
rm(data_to_remove)
rm(Data_nGy)


# 1-Build the Dose response to radiation without inhibitors
#####

#plot dose responses without inhibitors 
#(pool the 2 "0nM" points for BI907828 and for KTX169, 4N each = 16 points for each dose (Gy) for eaxh cell line)
#collect the 0nM points
response_all_0nM <- c()       # initialise aggregate of response_3_replicate for each cell line for the 2 ttt
#initialise iteration variables
i_ttt <- NULL
i_dose <- NULL
i_cell_line <- NULL
i_replicate <- NULL
for (i_ttt in c(1:length(ttt))) {   # ttt KTX-169 and BI907828
  for (i_dose in c(1:5)){           # doses Gy
    #initialise line number for each new dose (Gy) (i.e. each new plate)
    LNum <- 3 #first line used in 384-well plate
    for (i_cell_line in c(1:length(cell_lines))){ #cell lines
      #(re-)initialise response vectors
      response_4_replicate <- c()
      for (i_replicate in c(1:4)){ #replicates
        #fetch the line from raw data
        if (i_ttt==1) { #use either left or right half of the 384-well plate
          a <- 2
          b <- 12
        } else {
          a <- 13
          b <- 14
        }
        response_1_replicate <- c(unlist(RawData[[i_dose]][LNum, c(a, b)]))
        
        response_4_replicate <- c(response_4_replicate, response_1_replicate)
        #generate the name of the raw data line (n= 1 to n=3)
        name_response_4_replicate <- paste0(unlist(cell_lines[i_cell_line]), 
                                            "_", 
                                            unlist(ttt[i_ttt]), 
                                            "_0nM_",
                                            unlist(dose[i_dose]),
                                            "_4N")                   # eg. IB114_KTX169_0nM_0Gy_3N
        #assign data
        assign(x = name_response_4_replicate, value = response_4_replicate)
        #increment line number
        LNum <- LNum+1
        if (i_dose==1){                                                             # for the 0Gy doses, collect value for each cell line for standardisation
          name_cell_line_0nM_0Gy_STD <- paste0(unlist(cell_lines[i_cell_line]),
                                               "_",
                                               unlist(ttt[i_ttt]),
                                               "_0nM_",
                                               unlist(dose[i_dose]),
                                               "_STD")
          assign(x = name_cell_line_0nM_0Gy_STD, value = response_4_replicate)
        }
      }
      #name response vector with cell line name to retrieve for aggregation
      names(response_4_replicate) <- c(rep(unlist(cell_lines[i_cell_line]), 8))
      #aggregate all responses
      response_all_0nM <- c(response_all_0nM, response_4_replicate)
    }
  }
}


#build the matrix to plot

# pool 0nM from the 2 ttt
IB114_0nM_0Gy_STD <- c(IB114_KTX169_0nM_0Gy_STD, 
                      IB114_BI907828_0nM_0Gy_STD) 
IB115_0nM_0Gy_STD <- c(IB115_KTX169_0nM_0Gy_STD, 
                       IB115_BI907828_0nM_0Gy_STD)
IB128_0nM_0Gy_STD <- c(IB128_KTX169_0nM_0Gy_STD, 
                       IB128_BI907828_0nM_0Gy_STD)
list_0nM_0Gy_STD <- list(IB114_0nM_0Gy_STD, 
                         IB115_0nM_0Gy_STD, 
                         IB128_0nM_0Gy_STD)



for (i_cell_line in c(1:length(cell_lines))){
  cell_line_0nM <- subset(response_all_0nM, 
                          names(response_all_0nM)==unlist(cell_lines[i_cell_line])) # retrieve 0nM points for each cell line
  cell_line_0nM <- cell_line_0nM/mean(unlist(list_0nM_0Gy_STD[i_cell_line]))        # standardisation by the 0nM points at 0Gy
  
  #assign
  name_response_8_replicate <- paste0(unlist(cell_lines[i_cell_line]), "_0nM_AllDose")
  assign(x = name_response_8_replicate, value = cell_line_0nM)
}


dose_response_all_0nM <- rep(c(0, 2, 4, 6, 8), # 5 doses (Gy)
                             each = 8,         # 2 x 0nM points in quadriplicate = 8
                             times = 2)        # for the 2 ttt

data_0nM_AllDose <- as.data.frame(cbind(dose_response_all_0nM, 
                                        IB114_0nM_AllDose, 
                                        IB115_0nM_AllDose, 
                                        IB128_0nM_AllDose))
rownames(data_0nM_AllDose) <- NULL
data_0nM_AllDose <- data_0nM_AllDose[order(data_0nM_AllDose$dose_response_all_0nM), ]



#plot 0nM individually
ggplot(data = data_0nM_AllDose) +
  #scale_y_continuous(limits = c(0, 1.5)) +
  geom_point(aes(x = dose_response_all_0nM, y = IB114_0nM_AllDose),
             color = "darkgreen", size = 3, alpha = 0.3) +
  geom_smooth(aes(x = dose_response_all_0nM, y = IB114_0nM_AllDose),
              #method = "lm",
              color = "darkgreen") +
  scale_y_log10() 

ggplot(data = data_0nM_AllDose) +
  #scale_y_continuous(limits = c(0, 1.5)) +
  geom_point(aes(x = dose_response_all_0nM, y = IB115_0nM_AllDose),
             color = "red", size = 3, alpha = 0.3) +
  geom_smooth(aes(x = dose_response_all_0nM, y = IB115_0nM_AllDose),
              method = NULL,
              color = "red") +
  scale_y_log10() 

ggplot(data = data_0nM_AllDose) +
  geom_point(aes(x = dose_response_all_0nM, y = IB128_0nM_AllDose),
             color = "blue", size = 3, alpha = 0.3) +
  geom_smooth(aes(x = dose_response_all_0nM, y = IB128_0nM_AllDose),
              method = NULL,
              color = "blue") +
  scale_y_log10()


#plot 0nM for the 4 cell lines
ggplot(data = data_0nM_AllDose) +
  geom_point(aes(x = dose_response_all_0nM, y = IB114_0nM_AllDose), 
             color = "green", size = 3, alpha = 0.3) +
  geom_smooth(aes(x = dose_response_all_0nM, y = IB114_0nM_AllDose), 
              #method = "lm", 
              color = "green") +
  scale_y_log10() +
  
  geom_point(aes(x = dose_response_all_0nM, y = IB115_0nM_AllDose), 
             color = "red", size = 3, alpha = 0.3) +
  geom_smooth(aes(x = dose_response_all_0nM, y = IB115_0nM_AllDose), 
              #method = "lm", 
              color = "red") +
  scale_y_log10() +
  
  geom_point(aes(x = dose_response_all_0nM, y = IB128_0nM_AllDose), 
             color = "blue", size = 3, alpha = 0.3) +
  geom_smooth(aes(x = dose_response_all_0nM, y = IB128_0nM_AllDose),
              # method = "lm",
              color = "blue") +
  scale_y_log10() 



#make DRC dose-response curve of KTX169 and BI907828 for each dose

# Dose-response curve (nM) for each dose (Gy)
IC50s_ttt_dose <- c()         # Initialise IC50 table
IC80s_ttt_dose <- c()

# Dose-response curve (nM) for each dose (Gy)
# standardised by 0nM values of each dose Gy
# Dose-response curve (nM) for each dose (Gy) -----------------------------
# standardised by 0nM values of each dose Gy

#initialise iteration variables
i_ttt <- NULL
i_dose <- NULL
i_cell_line <- NULL
i_replicate <- NULL
#generate tables
for (i_ttt in c(1:length(ttt))) {            # for each of the 2 ttt
  #fetch the line from raw data
  if (i_ttt==1) {                      # use either left or right half of the 384-well plate
    Conc_range <- Conc_KTX169          # concentration ranges are different for each inhibitor
    a <- 2              # range of column to fetch data for KTX169
    b <- 12             # range of column to fetch data for KTX169
    c <- 2              # first 0nM point for for KTX169
    d <- 12             # second 0nM point for for KTX169
    matrix_ttt <- as.data.frame(Conc_KTX169, ncol = 1)       # initialise "big" matrix for each ttt
    matrix_ttt <- matrix_ttt[order(matrix_ttt[[1]]), , drop = FALSE]
    colnames(matrix_ttt) <- c("KTX169")
  } else {
    Conc_range <- Conc_BI907828
    a <- 13             # range of column to fetch data for BI907828
    b <- 23             # range of column to fetch data for BI907828
    c <- 13             # first 0nM point for for BI907828
    d <- 14             # second 0nM point for for BI907828
    matrix_ttt <- as.data.frame(Conc_BI907828, ncol = 1)       # initialise "big" matrix for each ttt
    matrix_ttt <- matrix_ttt[order(matrix_ttt[[1]]), , drop = FALSE]
    colnames(matrix_ttt) <- c("BI907828")
  }  # else
  
  rownames(matrix_ttt) <- NULL                        
  
  for (i_dose in c(1:length(dose))){                    # for each dose Gy
    #initialise line number for each new dose (i.e. each new plate)
    LNum <- 3                                  # first line used in 384-well plate
    
    
    for (i_cell_line in c(1:length(cell_lines))){             # for each cell line
      #re-initialise response vectors
      resp_col_4N <- c()                       # Initialise variable for 4 lines (4 replicates)
      average_0nM <- c()                       # Initialise mean 0nM ttt for one line (1 replicate)
      
      for (i_replicate in c(1:4)){           # replicates
        resp_col <- c(unlist(RawData[[i_dose]][LNum, c(a:b)]))  
        average_0nM <- mean(c(unlist(RawData[[i_dose]][LNum, c(c, d)])))            # normalised by each dose ttt=0nM
        resp_col <- resp_col/average_0nM                                            # normalised by each dose ttt=0nM
        
        
        
        #aggregate the 3n into one single vector. this will be the response column
        resp_col_4N <- c(resp_col_4N, resp_col)
        #assign data
        name_cell_lines_ttt_dose_4N <- paste0(unlist(cell_lines[i_cell_line]),      # generate the name of the raw data  for 4 lines
                                              "_",                                  # eg. IB111_KTX169_0Gy_4N
                                              unlist(ttt[i_ttt]),
                                              "_",
                                              unlist(dose[i_dose]),
                                              "_4N")
        assign(x = name_cell_lines_ttt_dose_4N, value = resp_col_4N)
        
        #iterate line number
        LNum <- LNum+1
      }  # replicate
      
      #Built data frame for each cell line
      matrix <- as.data.frame(Conc_range, nrow = 1)   # concentration ranges are different for each inhibitor
      matrix <- cbind(matrix, resp_col_4N)
      matrix <- matrix[order(matrix[[1]]), ]
      colnames(matrix) <- c(colnames(matrix)[1], 
                            paste0("col_", 
                                   name_cell_lines_ttt_dose_4N))
      rownames(matrix) <- NULL
      #assign matrix (re-initialised for each cell line)
      name_cell_lines_ttt_dose_4N_matrix <- paste0(name_cell_lines_ttt_dose_4N, 
                                                   "_matrix")
      
      assign(x = name_cell_lines_ttt_dose_4N_matrix, value = matrix)
      
      #built ttt matrix ("big" matrix) (re-initialised for each ttt)
      matrix_ttt <- cbind(matrix_ttt, matrix[[2]])
      colnames(matrix_ttt) <- c(colnames(matrix_ttt)[1:c(length(colnames(matrix_ttt))-1)], 
                                paste0("col_",
                                       name_cell_lines_ttt_dose_4N))
      
      
      
      #plot raw table check if works!!!!!!!does not seem to work!!!!!!!!!!!
      x <- unlist(matrix[1])             # x = resp_4N
      y <- unlist(matrix[2])             # y = concentration of ttt
      #ggplot(data = matrix) +
      #geom_point(aes(x = x, y = y)) +
      #scale_x_log10()
      
      #fit the model (4-parameters log-logistic model
      try(LL4_model <- drm(y ~ x, data = matrix,            # resp_3N ~ concentration of ttt
                           fct = LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50"),
                                      fixed = c(NA, NA, NA, NA))))
      #generate LL4_model name
      name_cell_lines_ttt_dose_LL4 <- paste0(unlist(cell_lines[i_cell_line]),      # generate the name of the raw data  for 4 lines
                                             "_",                                  # eg. IB111_KTX169_0Gy_4N
                                             unlist(ttt[i_ttt]),
                                             "_",
                                             unlist(dose[i_dose]),
                                             "_LL4")
      #assign LL4_model
      assign(x = name_cell_lines_ttt_dose_LL4, value = LL4_model)                  # assign LL4-model to plot it later
      #get IC50
      IC50 <- ED(LL4_model, c(50), interval = "delta")
      rownames(IC50) <- name_cell_lines_ttt_dose_4N
      IC50s_ttt_dose <- rbind(IC50s_ttt_dose, IC50)
      #get IC80
      IC80 <- ED(LL4_model, c(80), interval = "delta")
      rownames(IC80) <- name_cell_lines_ttt_dose_4N
      IC80s_ttt_dose <- rbind(IC80s_ttt_dose, IC80)
      #plot fitted model
      #plot(LL4_model, broken = TRUE, bp = 0.00001, type = "all", main = name_cell_lines_ttt_dose_3N, sub = IC50[1:1])
      #plot modeled response with confidence interval
      #confint(LL4_model)
      #plot(LL4_model, broken = TRUE, type = "confidence", main = name_cell_lines_ttt_dose_3N, sub = IC50[1:1])
      
      
    }   # cell line
  }     # dose Gy
  # assign "big" matrix_ttt
  name_matrix_ttt <- paste0("matrix_",
                            unlist(ttt[i_ttt]))
  
  assign(x = name_matrix_ttt, value = matrix_ttt)
  
}       # ttt
# end loop




# transform the 0nM in matrix into 1/100 of the 2nd lowest concentration
# To avoid: Warning message: In scale_x_log10() :   log-10 transformation introduced infinite values.
# the 0nM point can not be plotted in ggplot2 (div by 0, same issue as in graphpad prism)
matrix_KTX169_No_0 <- matrix_KTX169 %>% 
  mutate(KTX169 = 
           ifelse(matrix_KTX169$KTX169 == 0, 
                  yes = unique(matrix_KTX169$KTX169)[[which.min(sort(unique(matrix_KTX169$KTX169)))+1]]/100, 
                  no = matrix_KTX169$KTX169
           )
  )

matrix_BI907828_No_0 <- matrix_BI907828 %>% 
  mutate(BI907828 = 
           ifelse(matrix_BI907828$BI907828 == 0, 
                  yes = unique(matrix_BI907828$BI907828)[[which.min(sort(unique(matrix_BI907828$BI907828)))+1]]/100, 
                  no = matrix_BI907828$BI907828
           )
  )


# raw data from 2 surv curv
ggplot(data = matrix_BI907828_No_0) +
  geom_point(aes(x = BI907828, y = col_IB128_BI907828_0Gy_4N), color = "red") +
  geom_point(aes(x = BI907828, y = col_IB128_BI907828_2Gy_4N), color = "blue") +
  scale_x_log10()

# 1 data set: raw data + drm line
ggplot(data = matrix_KTX169_No_0, colour = "blue") +
  aes(x = KTX169, y = col_IB114_KTX169_0Gy_4N) + 
  geom_point(color = "blue") +
  geom_smooth(method = drm, method.args = list(fct = L.4()), se = FALSE) +
  scale_x_log10()

# KTX169
# se = FALSE ignores the confidence interval which is not handles by geom_smooth
ggplot(data = matrix_KTX169_No_0) +
  geom_smooth(aes(x = KTX169, y = col_IB114_KTX169_0Gy_4N, color = "0Gy"), colour = "black", method = drm, method.args = list(fct = L.4()), se = FALSE) +
  geom_point(aes(x = KTX169, y = col_IB114_KTX169_0Gy_4N), colour = "black", alpha = 0.2) +
  geom_smooth(aes(x = KTX169, y = col_IB114_KTX169_2Gy_4N, color = "2Gy"), colour = "green", method = drm, method.args = list(fct = L.4()), se = FALSE) +
  geom_point(aes(x = KTX169, y = col_IB114_KTX169_2Gy_4N), colour = "green", alpha = 0.2) +
  geom_smooth(aes(x = KTX169, y = col_IB114_KTX169_4Gy_4N, color = "4Gy"), colour = "blue", method = drm, method.args = list(fct = L.4()), se = FALSE) +
  geom_point(aes(x = KTX169, y = col_IB114_KTX169_4Gy_4N), colour = "blue", alpha = 0.2) +
  geom_smooth(aes(x = KTX169, y = col_IB114_KTX169_6Gy_4N, color = "6Gy"), colour = "grey", method = drm, method.args = list(fct = L.4()), se = FALSE) +
  geom_point(aes(x = KTX169, y = col_IB114_KTX169_6Gy_4N), colour = "grey", alpha = 0.2) +
  geom_smooth(aes(x = KTX169, y = col_IB114_KTX169_8Gy_4N, color = "8Gy"),colour = "red", method = drm, method.args = list(fct = L.4()), se = FALSE) +
  geom_point(aes(x = KTX169, y = col_IB114_KTX169_8Gy_4N), colour = "red", alpha = 0.2) +
  scale_x_log10()


ggplot(data = matrix_KTX169_No_0) +
  geom_smooth(aes(x = KTX169, y = col_IB115_KTX169_0Gy_4N, color = "0Gy"), colour = "black", method = drm, method.args = list(fct = L.4()), se = FALSE) +
  geom_point(aes(x = KTX169, y = col_IB115_KTX169_0Gy_4N), colour = "black", alpha = 0.2) +
  geom_smooth(aes(x = KTX169, y = col_IB115_KTX169_2Gy_4N, color = "2Gy"), colour = "green", method = drm, method.args = list(fct = L.4()), se = FALSE) +
  geom_point(aes(x = KTX169, y = col_IB115_KTX169_2Gy_4N), colour = "green", alpha = 0.2) +
  geom_smooth(aes(x = KTX169, y = col_IB115_KTX169_4Gy_4N, color = "4Gy"), colour = "blue", method = drm, method.args = list(fct = L.4()), se = FALSE) +
  geom_point(aes(x = KTX169, y = col_IB115_KTX169_4Gy_4N), colour = "blue", alpha = 0.2) +
  geom_smooth(aes(x = KTX169, y = col_IB115_KTX169_6Gy_4N, color = "6Gy"), colour = "grey", method = drm, method.args = list(fct = L.4()), se = FALSE) +
  geom_point(aes(x = KTX169, y = col_IB115_KTX169_6Gy_4N), colour = "grey", alpha = 0.2) +
  geom_smooth(aes(x = KTX169, y = col_IB115_KTX169_8Gy_4N, color = "8Gy"),colour = "red", method = drm, method.args = list(fct = L.4()), se = FALSE) +
  geom_point(aes(x = KTX169, y = col_IB115_KTX169_8Gy_4N), colour = "red", alpha = 0.2) +
  scale_x_log10()


ggplot(data = matrix_KTX169_No_0) +
  geom_smooth(aes(x = KTX169, y = col_IB128_KTX169_0Gy_4N, color = "0Gy"), colour = "black", method = drm, method.args = list(fct = L.4()), se = FALSE) +
  geom_point(aes(x = KTX169, y = col_IB128_KTX169_0Gy_4N), colour = "black", alpha = 0.2) +
  geom_smooth(aes(x = KTX169, y = col_IB128_KTX169_2Gy_4N, color = "2Gy"), colour = "green", method = drm, method.args = list(fct = L.4()), se = FALSE) +
  geom_point(aes(x = KTX169, y = col_IB128_KTX169_2Gy_4N), colour = "green", alpha = 0.2) +
  geom_smooth(aes(x = KTX169, y = col_IB128_KTX169_4Gy_4N, color = "4Gy"), colour = "blue", method = drm, method.args = list(fct = L.4()), se = FALSE) +
  geom_point(aes(x = KTX169, y = col_IB128_KTX169_4Gy_4N), colour = "blue", alpha = 0.2) +
  geom_smooth(aes(x = KTX169, y = col_IB128_KTX169_6Gy_4N, color = "6Gy"), colour = "grey", method = drm, method.args = list(fct = L.4()), se = FALSE) +
  geom_point(aes(x = KTX169, y = col_IB128_KTX169_6Gy_4N), colour = "grey", alpha = 0.2) +
  geom_smooth(aes(x = KTX169, y = col_IB128_KTX169_8Gy_4N, color = "8Gy"),colour = "red", method = drm, method.args = list(fct = L.4()), se = FALSE) +
  geom_point(aes(x = KTX169, y = col_IB128_KTX169_8Gy_4N), colour = "red", alpha = 0.2) +
  scale_x_log10()


# BI907828
ggplot(data = matrix_BI907828_No_0) +
  geom_smooth(aes(x = BI907828, y = col_IB114_BI907828_0Gy_4N, color = "0Gy"), colour = "black", method = drm, method.args = list(fct = L.4()), se = FALSE) +
  geom_point(aes(x = BI907828, y = col_IB114_BI907828_0Gy_4N), colour = "black", alpha = 0.2) +
  geom_smooth(aes(x = BI907828, y = col_IB114_BI907828_2Gy_4N, color = "2Gy"), colour = "green", method = drm, method.args = list(fct = L.4()), se = FALSE) +
  geom_point(aes(x = BI907828, y = col_IB114_BI907828_2Gy_4N), colour = "green", alpha = 0.2) +
  geom_smooth(aes(x = BI907828, y = col_IB114_BI907828_4Gy_4N, color = "4Gy"), colour = "blue", method = drm, method.args = list(fct = L.4()), se = FALSE) +
  geom_point(aes(x = BI907828, y = col_IB114_BI907828_4Gy_4N), colour = "blue", alpha = 0.2) +
  geom_smooth(aes(x = BI907828, y = col_IB114_BI907828_6Gy_4N, color = "6Gy"), colour = "grey", method = drm, method.args = list(fct = L.4()), se = FALSE) +
  geom_point(aes(x = BI907828, y = col_IB114_BI907828_6Gy_4N), colour = "grey", alpha = 0.2) +
  geom_smooth(aes(x = BI907828, y = col_IB114_BI907828_8Gy_4N, color = "8Gy"),colour = "red", method = drm, method.args = list(fct = L.4()), se = FALSE) +
  geom_point(aes(x = BI907828, y = col_IB114_BI907828_8Gy_4N), colour = "red", alpha = 0.2) +
  scale_x_log10()


ggplot(data = matrix_BI907828_No_0) +
  geom_smooth(aes(x = BI907828, y = col_IB115_BI907828_0Gy_4N, color = "0Gy"), colour = "black", method = drm, method.args = list(fct = L.4()), se = FALSE) +
  geom_point(aes(x = BI907828, y = col_IB115_BI907828_0Gy_4N), colour = "black", alpha = 0.2) +
  geom_smooth(aes(x = BI907828, y = col_IB115_BI907828_2Gy_4N, color = "2Gy"), colour = "green", method = drm, method.args = list(fct = L.4()), se = FALSE) +
  geom_point(aes(x = BI907828, y = col_IB115_BI907828_2Gy_4N), colour = "green", alpha = 0.2) +
  geom_smooth(aes(x = BI907828, y = col_IB115_BI907828_4Gy_4N, color = "4Gy"), colour = "blue", method = drm, method.args = list(fct = L.4()), se = FALSE) +
  geom_point(aes(x = BI907828, y = col_IB115_BI907828_4Gy_4N), colour = "blue", alpha = 0.2) +
  geom_smooth(aes(x = BI907828, y = col_IB115_BI907828_6Gy_4N, color = "6Gy"), colour = "grey", method = drm, method.args = list(fct = L.4()), se = FALSE) +
  geom_point(aes(x = BI907828, y = col_IB115_BI907828_6Gy_4N), colour = "grey", alpha = 0.2) +
  geom_smooth(aes(x = BI907828, y = col_IB115_BI907828_8Gy_4N, color = "8Gy"),colour = "red", method = drm, method.args = list(fct = L.4()), se = FALSE) +
  geom_point(aes(x = BI907828, y = col_IB115_BI907828_8Gy_4N), colour = "red", alpha = 0.2) +
  scale_x_log10()


ggplot(data = matrix_BI907828_No_0) +
  geom_smooth(aes(x = BI907828, y = col_IB128_BI907828_0Gy_4N, color = "0Gy"), colour = "black", method = drm, method.args = list(fct = L.4()), se = FALSE) +
  geom_point(aes(x = BI907828, y = col_IB128_BI907828_0Gy_4N), colour = "black", alpha = 0.2) +
  geom_smooth(aes(x = BI907828, y = col_IB128_BI907828_2Gy_4N, color = "2Gy"), colour = "green", method = drm, method.args = list(fct = L.4()), se = FALSE) +
  geom_point(aes(x = BI907828, y = col_IB128_BI907828_2Gy_4N), colour = "green", alpha = 0.2) +
  geom_smooth(aes(x = BI907828, y = col_IB128_BI907828_4Gy_4N, color = "4Gy"), colour = "blue", method = drm, method.args = list(fct = L.4()), se = FALSE) +
  geom_point(aes(x = BI907828, y = col_IB128_BI907828_4Gy_4N), colour = "blue", alpha = 0.2) +
  geom_smooth(aes(x = BI907828, y = col_IB128_BI907828_6Gy_4N, color = "6Gy"), colour = "grey", method = drm, method.args = list(fct = L.4()), se = FALSE) +
  geom_point(aes(x = BI907828, y = col_IB128_BI907828_6Gy_4N), colour = "grey", alpha = 0.2) +
  geom_smooth(aes(x = BI907828, y = col_IB128_BI907828_8Gy_4N, color = "8Gy"),colour = "red", method = drm, method.args = list(fct = L.4()), se = FALSE) +
  geom_point(aes(x = BI907828, y = col_IB128_BI907828_8Gy_4N), colour = "red", alpha = 0.2) +
  scale_x_log10()





########plot dose-response normalised with 0nM_0Gy
######
########plot dose-response normalised with 0nM_0Gy

#generate average 0nM 0Gy for each cell line (combine 0nM 0Gy for the 2 ttt (KTX and BI))
i_ttt <- NULL
i_dose <- NULL
i_cell_line <- NULL
i_replicate <- NULL
LNum <- 3
for (i_cell_line in c(1:length(cell_lines))) {
  points_0nM_0Gy <- c()
  for (i_replicate in c(1:replicates)) {
    points_0nM_0Gy_N <- unlist(RawData[["0Gy"]][LNum, c(2, 12, 13, 14)])
    points_0nM_0Gy <- c(points_0nM_0Gy, points_0nM_0Gy_N)
    LNum <- LNum + 1
  }
  name_points_0nM_0Gy <- paste0("data_0nM_0Gy_", unlist(cell_lines[i_cell_line]))
  assign(x = name_points_0nM_0Gy, value = points_0nM_0Gy)
}


#generate tables
#initialise iteration variables
i_ttt <- NULL
i_dose <- NULL
i_cell_line <- NULL
i_replicate <- NULL

for (i_ttt in c(1:length(ttt))) {            # for each of the 2 ttt
  #fetch the line from raw data
  if (i_ttt==1) {                      # use either left or right half of the 384-well plate
    Conc_range <- Conc_KTX169          # concentration ranges are different for each inhibitor
    a <- 2              # range of column to fetch data for KTX169
    b <- 12             # range of column to fetch data for KTX169
    c <- 2              # first 0nM point for for KTX169
    d <- 12             # second 0nM point for for KTX169
    matrix_ttt <- as.data.frame(Conc_KTX169, ncol = 1)       # initialise "big" matrix for each ttt
    matrix_ttt <- matrix_ttt[order(matrix_ttt[[1]]), , drop = FALSE]
    colnames(matrix_ttt) <- c("KTX169")
  } else {
    Conc_range <- Conc_BI907828
    a <- 13             # range of column to fetch data for BI907828
    b <- 23             # range of column to fetch data for BI907828
    c <- 13             # first 0nM point for for BI907828
    d <- 14             # second 0nM point for for BI907828
    matrix_ttt <- as.data.frame(Conc_BI907828, ncol = 1)       # initialise "big" matrix for each ttt
    matrix_ttt <- matrix_ttt[order(matrix_ttt[[1]]), , drop = FALSE]
    colnames(matrix_ttt) <- c("BI907828")
  }  # else
  
  rownames(matrix_ttt) <- NULL                        
  
  for (i_dose in c(1:length(dose))){                    # for each dose Gy
    #initialise line number for each new dose (i.e. each new plate)
    LNum <- 3                                  # first line used in 384-well plate
    
    for (i_cell_line in c(1:length(cell_lines))){             # for each cell line
      #re-initialise response vectors
      resp_col_4N_std0Gy <- c()                       # Initialise variable for 4 lines (4 replicates)
      
      
      for (i_replicate in c(1:replicates)){           # 4 replicates
        resp_col_std0Gy <- c(unlist(RawData[[i_dose]][LNum, c(a:b)]))  
        resp_col_std0Gy <- resp_col_std0Gy/mean(get(paste0("data_0nM_0Gy_", unlist(cell_lines[i_cell_line]))))*100              # normalised by ttt=0nM at 0Gy
        
        
        #aggregate the 3n into one single vector. this will be the response column
        resp_col_4N_std0Gy <- c(resp_col_4N_std0Gy, resp_col_std0Gy)
        #assign data
        name_cell_lines_ttt_dose_4N_std0Gy <- paste0(unlist(cell_lines[i_cell_line]),      # generate the name of the raw data  for 4 lines
                                                     "_",                                  # eg. IB111_KTX169_0Gy_4N
                                                     unlist(ttt[i_ttt]),
                                                     "_",
                                                     unlist(dose[i_dose]),
                                                     "_4N_std0Gy")
        assign(x = name_cell_lines_ttt_dose_4N_std0Gy, value = resp_col_4N_std0Gy)
        
        #iterate line number
        LNum <- LNum+1
      }  # replicate
      
      #Built data frame for each cell line
      matrix <- as.data.frame(Conc_range, nrow = 1)   # concentration ranges are different for each inhibitor
      matrix <- cbind(matrix, resp_col_4N_std0Gy)
      matrix <- matrix[order(matrix[[1]]), ]
      colnames(matrix) <- c(colnames(matrix)[1], 
                            paste0("col_", 
                                   name_cell_lines_ttt_dose_4N_std0Gy))
      rownames(matrix) <- NULL
      #assign matrix
      name_cell_lines_ttt_dose_4N_matrix_std0Gy <- paste0(name_cell_lines_ttt_dose_4N_std0Gy, 
                                                          "_matrix")
      
      assign(x = name_cell_lines_ttt_dose_4N_matrix_std0Gy, value = matrix)
      
      #built ttt matrix ("big" matrix)
      matrix_ttt <- cbind(matrix_ttt, matrix[[2]])
      colnames(matrix_ttt) <- c(colnames(matrix_ttt)[1:c(length(colnames(matrix_ttt))-1)], 
                                paste0("std0Gy_",
                                       name_cell_lines_ttt_dose_4N_std0Gy))
      
      
      
      #plot raw table check if works!!!!!!!does not seem to work!!!!!!!!!!!
      x <- unlist(matrix[1])             # x = resp_4N
      y <- unlist(matrix[2])             # y = concentration of ttt
      #ggplot(data = matrix) +
      #geom_point(aes(x = x, y = y)) +
      #scale_x_log10()
      
      #fit the model (4-parameters log-logistic model
      try(LL4_model <- drm(y ~ x, data = matrix,            # resp_4N ~ concentration of ttt
                           fct = LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50"),
                                      fixed = c(NA, NA, NA, NA))))
      #generate LL4_model name
      name_cell_lines_ttt_dose_LL4_std0Gy <- paste0(unlist(cell_lines[i_cell_line]),      # generate the name of the raw data  for 4 lines
                                                    "_",                                  # eg. IB111_KTX169_0Gy_4N
                                                    unlist(ttt[i_ttt]),
                                                    "_",
                                                    unlist(dose[i_dose]),
                                                    "_LL4_std0Gy")
      #assign LL4_model
      assign(x = name_cell_lines_ttt_dose_LL4_std0Gy, value = LL4_model)                  # assign LL4-model to plot it later
      #get IC50
      IC50 <- ED(LL4_model, c(50), interval = "delta")
      rownames(IC50) <- name_cell_lines_ttt_dose_4N_std0Gy
      IC50s_ttt_dose_std0Gy <- rbind(IC50s_ttt_dose, IC50)
      #get IC80
      IC80 <- ED(LL4_model, c(80), interval = "delta")
      rownames(IC80) <- name_cell_lines_ttt_dose_4N_std0Gy
      IC80s_ttt_dose_std0Gy <- rbind(IC80s_ttt_dose, IC80)
      #plot fitted model
      #plot(LL4_model, broken = TRUE, bp = 0.00001, type = "all", main = name_cell_lines_ttt_dose_3N, sub = IC50[1:1])
      #plot modeled response with confidence interval
      #confint(LL4_model)
      #plot(LL4_model, broken = TRUE, type = "confidence", main = name_cell_lines_ttt_dose_3N, sub = IC50[1:1])
      
      
    }   # cell line
  }     # dose Gy
  # assign "big" matrix_ttt
  name_matrix_ttt_std0Gy <- paste0("matrix_",
                                   unlist(ttt[i_ttt]),
                                   "_std0Gy")
  
  assign(x = name_matrix_ttt_std0Gy, value = matrix_ttt)
  
}       # ttt
# end loop




# transform the 0nM in matrix into 1/100 of the 2nd lowest concentration
# To avoid: Warning message: In scale_x_log10() :   log-10 transformation introduced infinite values.
# the 0nM point can not be plotted in ggplot2 (div by 0, same issue as in graphpad prism)
matrix_KTX169_No_0_std0Gy <- matrix_KTX169_std0Gy %>% 
  mutate(KTX169 = 
           ifelse(matrix_KTX169_std0Gy$KTX169 == 0, 
                  yes = unique(matrix_KTX169_std0Gy$KTX169)[[which.min(sort(unique(matrix_KTX169_std0Gy$KTX169)))+1]]/100, 
                  no = matrix_KTX169_std0Gy$KTX169
           )
  )

matrix_BI907828_No_0_std0Gy <- matrix_BI907828_std0Gy %>% 
  mutate(BI907828 = 
           ifelse(matrix_BI907828_std0Gy$BI907828 == 0, 
                  yes = unique(matrix_BI907828_std0Gy$BI907828)[[which.min(sort(unique(matrix_BI907828_std0Gy$BI907828)))+1]]/100, 
                  no = matrix_BI907828_std0Gy$BI907828
           )
  )


# raw data from 2 surv curv
ggplot(data = matrix_KTX169_No_0_std0Gy) +
  geom_point(aes(x = KTX169, y = std0Gy_IB114_KTX169_0Gy_4N_std0Gy), color = "red") +
  geom_point(aes(x = KTX169, y = std0Gy_IB114_KTX169_2Gy_4N_std0Gy), color = "blue") +
  scale_x_log10()

# 1 data set: raw data + drm line
ggplot(data = matrix_KTX169_No_0_std0Gy, colour = "blue") +
  aes(x = KTX169, y = std0Gy_IB114_KTX169_0Gy_4N_std0Gy) + 
  geom_point(color = "blue") +
  geom_smooth(method = drm, method.args = list(fct = L.4()), se = FALSE) +
  scale_x_log10()

# KTX169
ggplot(data = matrix_KTX169_No_0_std0Gy) +
  geom_smooth(aes(x = KTX169, y = std0Gy_IB114_KTX169_0Gy_4N_std0Gy, color = "0Gy"), colour = "black", method = drm, method.args = list(fct = L.4()), se = FALSE) +
  geom_point(aes(x = KTX169, y = std0Gy_IB114_KTX169_0Gy_4N_std0Gy), colour = "black", alpha = 0.2) +
  geom_smooth(aes(x = KTX169, y = std0Gy_IB114_KTX169_2Gy_4N_std0Gy, color = "2Gy"), colour = "green", method = drm, method.args = list(fct = L.4()), se = FALSE) +
  geom_point(aes(x = KTX169, y = std0Gy_IB114_KTX169_2Gy_4N_std0Gy), colour = "green", alpha = 0.2) +
  geom_smooth(aes(x = KTX169, y = std0Gy_IB114_KTX169_4Gy_4N_std0Gy, color = "4Gy"), colour = "blue", method = drm, method.args = list(fct = L.4()), se = FALSE) +
  geom_point(aes(x = KTX169, y = std0Gy_IB114_KTX169_4Gy_4N_std0Gy), colour = "blue", alpha = 0.2) +
  geom_smooth(aes(x = KTX169, y = std0Gy_IB114_KTX169_6Gy_4N_std0Gy, color = "6Gy"), colour = "grey", method = drm, method.args = list(fct = L.4()), se = FALSE) +
  geom_point(aes(x = KTX169, y = std0Gy_IB114_KTX169_6Gy_4N_std0Gy), colour = "grey", alpha = 0.2) +
  geom_smooth(aes(x = KTX169, y = std0Gy_IB114_KTX169_8Gy_4N_std0Gy, color = "8Gy"),colour = "red", method = drm, method.args = list(fct = L.4()), se = FALSE) +
  geom_point(aes(x = KTX169, y = std0Gy_IB114_KTX169_8Gy_4N_std0Gy), colour = "red", alpha = 0.2) +
  scale_x_log10()


ggplot(data = matrix_KTX169_No_0_std0Gy) +
  geom_smooth(aes(x = KTX169, y = std0Gy_IB115_KTX169_0Gy_4N_std0Gy, color = "0Gy"), colour = "black", method = drm, method.args = list(fct = L.4()), se = FALSE) +
  geom_point(aes(x = KTX169, y = std0Gy_IB115_KTX169_0Gy_4N_std0Gy), colour = "black", alpha = 0.2) +
  geom_smooth(aes(x = KTX169, y = std0Gy_IB115_KTX169_2Gy_4N_std0Gy, color = "2Gy"), colour = "green", method = drm, method.args = list(fct = L.4()), se = FALSE) +
  geom_point(aes(x = KTX169, y = std0Gy_IB115_KTX169_2Gy_4N_std0Gy), colour = "green", alpha = 0.2) +
  geom_smooth(aes(x = KTX169, y = std0Gy_IB115_KTX169_4Gy_4N_std0Gy, color = "4Gy"), colour = "blue", method = drm, method.args = list(fct = L.4()), se = FALSE) +
  geom_point(aes(x = KTX169, y = std0Gy_IB115_KTX169_4Gy_4N_std0Gy), colour = "blue", alpha = 0.2) +
  geom_smooth(aes(x = KTX169, y = std0Gy_IB115_KTX169_6Gy_4N_std0Gy, color = "6Gy"), colour = "grey", method = drm, method.args = list(fct = L.4()), se = FALSE) +
  geom_point(aes(x = KTX169, y = std0Gy_IB115_KTX169_6Gy_4N_std0Gy), colour = "grey", alpha = 0.2) +
  geom_smooth(aes(x = KTX169, y = std0Gy_IB115_KTX169_8Gy_4N_std0Gy, color = "8Gy"),colour = "red", method = drm, method.args = list(fct = L.4()), se = FALSE) +
  geom_point(aes(x = KTX169, y = std0Gy_IB115_KTX169_8Gy_4N_std0Gy), colour = "red", alpha = 0.2) +
  scale_x_log10()



ggplot(data = matrix_KTX169_No_0_std0Gy) +
  geom_smooth(aes(x = KTX169, y = std0Gy_IB128_KTX169_0Gy_4N_std0Gy, color = "0Gy"), colour = "black", method = drm, method.args = list(fct = L.4()), se = FALSE) +
  geom_point(aes(x = KTX169, y = std0Gy_IB128_KTX169_0Gy_4N_std0Gy), colour = "black", alpha = 0.2) +
  geom_smooth(aes(x = KTX169, y = std0Gy_IB128_KTX169_2Gy_4N_std0Gy, color = "2Gy"), colour = "green", method = drm, method.args = list(fct = L.4()), se = FALSE) +
  geom_point(aes(x = KTX169, y = std0Gy_IB128_KTX169_2Gy_4N_std0Gy), colour = "green", alpha = 0.2) +
  geom_smooth(aes(x = KTX169, y = std0Gy_IB128_KTX169_4Gy_4N_std0Gy, color = "4Gy"), colour = "blue", method = drm, method.args = list(fct = L.4()), se = FALSE) +
  geom_point(aes(x = KTX169, y = std0Gy_IB128_KTX169_4Gy_4N_std0Gy), colour = "blue", alpha = 0.2) +
  geom_smooth(aes(x = KTX169, y = std0Gy_IB128_KTX169_6Gy_4N_std0Gy, color = "6Gy"), colour = "grey", method = drm, method.args = list(fct = L.4()), se = FALSE) +
  geom_point(aes(x = KTX169, y = std0Gy_IB128_KTX169_6Gy_4N_std0Gy), colour = "grey", alpha = 0.2) +
  geom_smooth(aes(x = KTX169, y = std0Gy_IB128_KTX169_8Gy_4N_std0Gy, color = "8Gy"),colour = "red", method = drm, method.args = list(fct = L.4()), se = FALSE) +
  geom_point(aes(x = KTX169, y = std0Gy_IB128_KTX169_8Gy_4N_std0Gy), colour = "red", alpha = 0.2) +
  scale_x_log10()


# BI907828
ggplot(data = matrix_BI907828_No_0_std0Gy) +
  geom_smooth(aes(x = BI907828, y = std0Gy_IB114_BI907828_0Gy_4N_std0Gy, color = "0Gy"), colour = "black", method = drm, method.args = list(fct = L.4()), se = FALSE) +
  geom_point(aes(x = BI907828, y = std0Gy_IB114_BI907828_0Gy_4N_std0Gy), colour = "black", alpha = 0.2) +
  geom_smooth(aes(x = BI907828, y = std0Gy_IB114_BI907828_2Gy_4N_std0Gy, color = "2Gy"), colour = "green", method = drm, method.args = list(fct = L.4()), se = FALSE) +
  geom_point(aes(x = BI907828, y = std0Gy_IB114_BI907828_2Gy_4N_std0Gy), colour = "green", alpha = 0.2) +
  geom_smooth(aes(x = BI907828, y = std0Gy_IB114_BI907828_4Gy_4N_std0Gy, color = "4Gy"), colour = "blue", method = drm, method.args = list(fct = L.4()), se = FALSE) +
  geom_point(aes(x = BI907828, y = std0Gy_IB114_BI907828_4Gy_4N_std0Gy), colour = "blue", alpha = 0.2) +
  geom_smooth(aes(x = BI907828, y = std0Gy_IB114_BI907828_6Gy_4N_std0Gy, color = "6Gy"), colour = "grey", method = drm, method.args = list(fct = L.4()), se = FALSE) +
  geom_point(aes(x = BI907828, y = std0Gy_IB114_BI907828_6Gy_4N_std0Gy), colour = "grey", alpha = 0.2) +
  geom_smooth(aes(x = BI907828, y = std0Gy_IB114_BI907828_8Gy_4N_std0Gy, color = "8Gy"),colour = "red", method = drm, method.args = list(fct = L.4()), se = FALSE) +
  geom_point(aes(x = BI907828, y = std0Gy_IB114_BI907828_8Gy_4N_std0Gy), colour = "red", alpha = 0.2) +
  scale_x_log10()


ggplot(data = matrix_BI907828_No_0_std0Gy) +
  geom_smooth(aes(x = BI907828, y = std0Gy_IB115_BI907828_0Gy_4N_std0Gy, color = "0Gy"), colour = "black", method = drm, method.args = list(fct = L.4()), se = FALSE) +
  geom_point(aes(x = BI907828, y = std0Gy_IB115_BI907828_0Gy_4N_std0Gy), colour = "black", alpha = 0.2) +
  geom_smooth(aes(x = BI907828, y = std0Gy_IB115_BI907828_2Gy_4N_std0Gy, color = "2Gy"), colour = "green", method = drm, method.args = list(fct = L.4()), se = FALSE) +
  geom_point(aes(x = BI907828, y = std0Gy_IB115_BI907828_2Gy_4N_std0Gy), colour = "green", alpha = 0.2) +
  geom_smooth(aes(x = BI907828, y = std0Gy_IB115_BI907828_4Gy_4N_std0Gy, color = "4Gy"), colour = "blue", method = drm, method.args = list(fct = L.4()), se = FALSE) +
  geom_point(aes(x = BI907828, y = std0Gy_IB115_BI907828_4Gy_4N_std0Gy), colour = "blue", alpha = 0.2) +
  geom_smooth(aes(x = BI907828, y = std0Gy_IB115_BI907828_6Gy_4N_std0Gy, color = "6Gy"), colour = "grey", method = drm, method.args = list(fct = L.4()), se = FALSE) +
  geom_point(aes(x = BI907828, y = std0Gy_IB115_BI907828_6Gy_4N_std0Gy), colour = "grey", alpha = 0.2) +
  geom_smooth(aes(x = BI907828, y = std0Gy_IB115_BI907828_8Gy_4N_std0Gy, color = "8Gy"),colour = "red", method = drm, method.args = list(fct = L.4()), se = FALSE) +
  geom_point(aes(x = BI907828, y = std0Gy_IB115_BI907828_8Gy_4N_std0Gy), colour = "red", alpha = 0.2) +
  scale_x_log10()



ggplot(data = matrix_BI907828_No_0_std0Gy) +
  geom_smooth(aes(x = BI907828, y = std0Gy_IB128_BI907828_0Gy_4N_std0Gy, color = "0Gy"), colour = "black", method = drm, method.args = list(fct = L.4()), se = FALSE) +
  geom_point(aes(x = BI907828, y = std0Gy_IB128_BI907828_0Gy_4N_std0Gy), colour = "black", alpha = 0.2) +
  geom_smooth(aes(x = BI907828, y = std0Gy_IB128_BI907828_2Gy_4N_std0Gy, color = "2Gy"), colour = "green", method = drm, method.args = list(fct = L.4()), se = FALSE) +
  geom_point(aes(x = BI907828, y = std0Gy_IB128_BI907828_2Gy_4N_std0Gy), colour = "green", alpha = 0.2) +
  geom_smooth(aes(x = BI907828, y = std0Gy_IB128_BI907828_4Gy_4N_std0Gy, color = "4Gy"), colour = "blue", method = drm, method.args = list(fct = L.4()), se = FALSE) +
  geom_point(aes(x = BI907828, y = std0Gy_IB128_BI907828_4Gy_4N_std0Gy), colour = "blue", alpha = 0.2) +
  geom_smooth(aes(x = BI907828, y = std0Gy_IB128_BI907828_6Gy_4N_std0Gy, color = "6Gy"), colour = "grey", method = drm, method.args = list(fct = L.4()), se = FALSE) +
  geom_point(aes(x = BI907828, y = std0Gy_IB128_BI907828_6Gy_4N_std0Gy), colour = "grey", alpha = 0.2) +
  geom_smooth(aes(x = BI907828, y = std0Gy_IB128_BI907828_8Gy_4N_std0Gy, color = "8Gy"),colour = "red", method = drm, method.args = list(fct = L.4()), se = FALSE) +
  geom_point(aes(x = BI907828, y = std0Gy_IB128_BI907828_8Gy_4N_std0Gy), colour = "red", alpha = 0.2) +
  scale_x_log10()






#class(IC50s_ttt_dose)
#IC80_to_print <- as.data.frame(IC80s_ttt_dose)
#rownames(IC80_to_print) 
#IC80_to_print <- cbind(rownames(IC80_to_print), IC80_to_print)
#write.table(x = IC80_to_print, file = "IC80.xls", sep = "\t")

#KTX <- grep(pattern = "KTX", x = IC80_to_print$`rownames(IC80_to_print)`)

