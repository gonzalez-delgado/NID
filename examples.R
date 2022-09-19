

# Path where the data for each central amino-acid is located
data_path <- getwd() 

# Amino-acid considered
amino_acid_list<-c("ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU","HIS", "ILE", "LEU", "LYS", "MET", "PHE", "SER", "THR", "TRP", "TYR", "VAL")

# Independence test (L and R independent for a given C and (phi, psi))

test_list_I <- independence_test(central_name = 'ALA', aa_list = amino_acid_list,
                   N_bin = 30, N_n = 200, ball = FALSE, plot = FALSE, N_rep = NULL, set_path = data_path) # Discretization method I

test_list_I[!is.na(test_list_I$p_value),] # Show results with admitted grid divisions


test_list_II <- independence_test(central_name = 'ALA', aa_list = amino_acid_list,
                     N_bin = 30, N_n = 1000, ball = TRUE, plot = FALSE, N_rep = NULL, dist_tol = 0.15, set_path = data_path) # Discretization method II

test_list_II[!is.na(test_list_II$p_value),] # Show results with admitted grid divisions


test_list_III <- independence_test(central_name = 'ALA', aa_list = amino_acid_list,
                      N_bin = 30, N_n = 1000, ball = TRUE, plot = FALSE, N_rep = 1000, dist_tol = 0.15, set_path = data_path) # Discretization method III

test_list_III[!is.na(test_list_III$p_value),] # Show results with admitted grid divisions


# Effect of polarity

amino_acid_list_polar <- c('ARG', 'LYS', 'ASP', 'GLU', 'ASN', 'GLN')
amino_acid_list_nonpolar <- c('ALA', 'ILE', 'LEU', 'MET', 'PHE', 'VAL')

polar_polar_test <- try(independence_test('ALA', aa_list = amino_acid_list, aa_list_left = amino_acid_list_polar, aa_list_right = amino_acid_list_polar, 
                                                          N_bin = 30, N_n = 1000, ball = TRUE, plot = FALSE, N_rep = 1000, dist_tol = 0.15, set_path = data_path))
  
nonpolar_nonpolar_test <- try(independence_test('ALA', aa_list = amino_acid_list, aa_list_left = amino_acid_list_nonpolar, aa_list_right = amino_acid_list_nonpolar,
                                                             N_bin = 30, N_n = 1000, ball = TRUE, plot = FALSE, N_rep = 100, dist_tol = 0.15, set_path = data_path))
  
polar_nonpolar_test <- try(independence_test('ALA', aa_list = amino_acid_list, aa_list_left = amino_acid_list_polar, aa_list_right = amino_acid_list_nonpolar,
                                                             N_bin = 30, N_n = 1000, ball = TRUE, plot = FALSE, N_rep = 1000, dist_tol = 0.15, set_path = data_path))
  
nonpolar_polar_test <- try(independence_test('ALA', aa_list = amino_acid_list, aa_list_left = amino_acid_list_nonpolar, aa_list_right = amino_acid_list_polar,
                                                             N_bin = 30, N_n = 1000, ball = TRUE, plot = FALSE, N_rep = 1000, dist_tol = 0.15, set_path = data_path))

# AUC computation

polar <- polar_polar_test[which(polar_polar_test$mat_dim == 6), ]
non_polar <- nonpolar_nonpolar_test[which(nonpolar_nonpolar_test$mat_dim == 6), ]
polar_nonpolar <- polar_nonpolar_test[which(polar_nonpolar_test$mat_dim == 6), ]
nonpolar_polar <- nonpolar_polar_test[which(nonpolar_polar_test$mat_dim == 6), ]
  
f_polar <- ecdf(polar$p_value)
f_nonpolar <- ecdf(non_polar$p_value)
f_polar_nonpolar <- ecdf(polar_nonpolar$p_value)
f_nonpolar_polar <- ecdf(nonpolar_polar$p_value)
  
auc_polar <- integrate(f_polar, 0, 1, subdivisions = 10^3)$value
auc_nonpolar <- integrate(f_nonpolar, 0, 1, subdivisions = 10^3)$value
auc_polar_nonpolar <- integrate(f_polar_nonpolar, 0, 1, subdivisions = 10^3)$value
auc_nonpolar_polar <- integrate(f_nonpolar_polar, 0, 1, subdivisions = 10^3)$value

# Show results
cat(paste('AUC (P-P): ',auc_polar, '\nAUC (H-H): ',auc_nonpolar,'\nAUC (P-H): ',auc_polar_nonpolar,'\nAUC (H-P): ',auc_nonpolar_polar))

