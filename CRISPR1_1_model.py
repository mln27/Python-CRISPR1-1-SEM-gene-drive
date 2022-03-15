# This module runs the CRISPR one-target location to one-target location using the equations generated and parameters
# set in the input file
#
# Model developed by Martial Ndeffo (m.ndeffo@tamu.edu)

import pandas as pd
import numpy as np
import math
import time
import CRISPR1_1_Equations
import matplotlib.pyplot as plt


def CRISPR_1_1_Progression():
    
    # stage duration and growth rate parameters
    T_e = 5
    T_l = 6
    T_p = 4
    r_m = 1.175
    mu_e = 0.175
    mu_p = 0.175

    cols = generate_allele_combinations(['w', 'v', 'u', 'r', 'g', 's'])
    results_temp = []

    # Reading and assigning values from Excel document
    M = pd.read_excel('CRISPR1-1 Input Parameters.xlsx')
    
    g0 = np.array(M.iloc[20:62, 3].values)

    sigma = M.iloc[5, 3]
    lamda = M.iloc[6, 3]

    q2 = M.iloc[8, 3]
    q1 = M.iloc[9, 3]
    P2 = M.iloc[10, 3]
    P1 = M.iloc[11, 3]
    delta2 = M.iloc[12, 3]
    delta1 = M.iloc[13, 3]
    mortRateAdult1 = np.array([M.iloc[14, 3]] * 2 * len(cols))
    mu_ad = mortRateAdult1[0]
    epsilon1 = M.iloc[16, 3]
    epsilon2 = M.iloc[15, 3]
    R_m = r_m**(T_e + T_l + T_p + 1/mu_ad)
    mu_l = 1 - (R_m*mu_ad/(0.5*lamda*(1-mu_ad)))**(1/(T_e + T_l + T_p))
    theta_e = (1-mu_e)**T_e
    theta_l = (1-mu_l)**T_l
    theta_p = (1-mu_p)**T_p
    alpha = (0.5*lamda*theta_e*sum(g0)/(R_m-1))*((1-theta_l/R_m)/(1-(theta_l/R_m)**(1/T_l)))
    mortRateJuven = np.array([mu_l] * 2 * len(cols))
    developTime = T_e + T_l + T_p
    gens = M.iloc[17, 3]
    span = (gens)*developTime
    

    male_alphas_init = M.iloc[0, 3:-1].values
    male_alphas = [i for i in male_alphas_init if not math.isnan(i)]  # refers to the male alpha values entered
    male_gammas_init = M.iloc[1, 3:-1].values
    male_gammas = [i for i in male_gammas_init if not math.isnan(i)]  # refers to the male gamma values entered

    female_alphas_init = M.iloc[2, 3:-1].values
    female_alphas = [i for i in female_alphas_init if not math.isnan(i)]  # refers to the male alpha values entered
    female_gammas_init = M.iloc[3, 3:-1].values
    female_gammas = [i for i in female_gammas_init if not math.isnan(i)]  # refers to the male gamma values entered

    # Assignment of initial fitness costs
    fitCost = np.array(M.iloc[20:62, 9].values)

    # recalculate fitness costs
    for index,val in enumerate(fitCost):
        mortRateAdult1[index] = mortRateAdult1[index] * (1 + fitCost[index])
     
    # calculate the time period which to track carry-out adults
    min_adult_mort1 = min(mortRateAdult1)
    
    
    # RUN THE MODEL
    pop = [0 for x in range(span)]

    # indexes are day
    surviving_Adults = [[0 for x in range(2 * len(cols))] for y in range(span)]
    surviving_Larvae = [[0 for x in range(2 * len(cols))] for y in range(span)]
    seeded_Larvae = [[0 for x in range(2 * len(cols))] for y in range(span)]
    seeded_Eggs = [[0 for x in range(2 * len(cols))] for y in range(span)]
    total_Adults = [[0 for x in range(2 * len(cols))] for y in range(span)]
    total_Larvae = [[0 for x in range(2 * len(cols))] for y in range(span)]
    total_Eggs = [[0 for x in range(2 * len(cols))] for y in range(span)]
    total_population_by_genotype = [[0 for x in range(2 * len(cols))] for y in range(span)]

    for alpha1 in male_alphas:
        for alpha2 in female_alphas:
            for gamma1 in male_gammas:
                for gamma2 in female_gammas:

                    beta1 = 1 - (alpha1 + gamma1)
                    beta2 = 1 - (alpha2 + gamma2)

                    # loop through each time step
                    start_time = time.time()
                    for T in range(span):
                        if T == 0:
                            total_Adults[0] = g0
                            total_Larvae[0][0] = alpha*(R_m - 1)/2
                            total_Larvae[0][21] = alpha*(R_m - 1)/2
                            proportionPop = convert_to_proportion(g0)
                            seeded_Eggs[0] = CRISPR1_1_Equations.CRISPR1_1_Rates(proportionPop, sigma, lamda, delta1, delta2, fitCost,
                                                                                      q1, q2, P1, P2, alpha1, alpha2, beta1, beta2, gamma1, gamma2, epsilon1, epsilon2)
                            total_Eggs[T] = [i + j for i, j in zip(total_Eggs[T], seeded_Eggs[T])]
                            
                            #Prior_Eggs is 0 if the population is assumed to be seeded without prior history
                            #rather than having naturally evolved to its equilibrium 
                            #g01=g0
                            #g01[18]=0
                            #proportionPop = convert_to_proportion(g01)
                            Prior_Eggs = [0 for x in range(2 * len(cols))]
                            #CRISPR1_1_Equations.CRISPR1_1_Rates(proportionPop, sigma, lamda, delta1, delta2, fitCost,
                            #                                                          q1, q2, P1, P2, alpha1, alpha2, beta1, beta2, gamma1, gamma2, epsilon1, epsilon2)
                        elif 0 < T < (T_e+1):

                            # only the initial adults will be producing offspring at this time.
                            surviving_Adults[T] = [i * j for i, j in zip(total_Adults[T-1], [(1 - k) for k in mortRateAdult1])]
                           
                            proportionPop = convert_to_proportion(surviving_Adults[T])
                            J2 = CRISPR1_1_Equations.CRISPR1_1_Rates(proportionPop, sigma, lamda, delta1, delta2, fitCost, q1, q2,
                                                                     P1, P2, alpha1, alpha2, beta1, beta2, gamma1, gamma2, epsilon1, epsilon2)
                            seeded_Eggs[T] = [i + j for i, j in zip(seeded_Eggs[T], J2)]
                            
                            F_prod = 1
                            for i in range(1,T_l+1):
                                F_prod = F_prod*F(alpha,sum(total_Larvae[0]),T_l)
                                
                            new_pupae = [i * j for i, j in zip(Prior_Eggs, [theta_e*theta_l*F_prod for k in range(len(Prior_Eggs))])]
                            surviving_Larvae[T] = [i*j for i, j in zip(total_Larvae[T-1], [(1-mu_l)*F(alpha,sum(total_Larvae[T-1]),T_l) for k in range(len(total_Larvae[T-1]))])]
                            seeded_Larvae[T] = [i * j for i, j in zip(Prior_Eggs, [theta_e for k in range(len(Prior_Eggs))])]
                            total_Larvae[T] = [i + j + k - q for i, j, k, q in zip(total_Larvae[T], surviving_Larvae[T], seeded_Larvae[T], new_pupae)]
                            # new adults create genotypes of
                            F_prod = 1
                            for i in range(1,T_l+1):
                                F_prod = F_prod*F(alpha,sum(total_Larvae[0]),T_l)
                                
                            new_adults = [i * j for i, j in zip(Prior_Eggs, [theta_e*theta_l*theta_p*F_prod*(1-mu_ad) for k in range(len(Prior_Eggs))])]
                            total_Adults[T] = [i + j + k for i, j, k in zip(total_Adults[T], surviving_Adults[T], new_adults)]
                            total_Eggs[T] = [i + j for i, j in zip(total_Eggs[T], J2)]
                           
                        
                        elif T_e < T < (T_e + T_l+1):

                             # only the initial adults will be producing offspring at this time.
                            surviving_Adults[T] = [i * j for i, j in zip(total_Adults[T-1], [(1 - k) for k in mortRateAdult1])]
                          
                            proportionPop = convert_to_proportion(surviving_Adults[T])
                            J2 = CRISPR1_1_Equations.CRISPR1_1_Rates(proportionPop, sigma, lamda, delta1, delta2, fitCost, q1, q2,
                                                                     P1, P2, alpha1, alpha2, beta1, beta2, gamma1, gamma2, epsilon1, epsilon2)
                            seeded_Eggs[T] = [i + j for i, j in zip(seeded_Eggs[T], J2)]
                            
                            F_prod = 1
                            for i in range(1,T_l+1):
                                F_prod = F_prod*F(alpha,sum(total_Larvae[0]),T_l)
                                
                            new_pupae = [i * j for i, j in zip(Prior_Eggs, [theta_e*theta_l*F_prod for k in range(len(Prior_Eggs))])]
                            surviving_Larvae[T] = [i*j for i, j in zip(total_Larvae[T-1], [(1-mu_l)*F(alpha,sum(total_Larvae[T-1]),T_l) for k in range(len(total_Larvae[T-1]))])]
                            seeded_Larvae[T] = [i * j for i, j in zip(seeded_Eggs[T - T_e], [theta_e for k in range(len(seeded_Eggs[T - T_e]))])]
                            total_Larvae[T] = [i + j + k - q for i, j, k, q in zip(total_Larvae[T], surviving_Larvae[T], seeded_Larvae[T], new_pupae)]
                            # new adults create genotypes of
                            F_prod = 1
                            for i in range(1,T_l+1):
                                F_prod = F_prod*F(alpha,sum(total_Larvae[0]),T_l)
                                
                            new_adults = [i * j for i, j in zip(Prior_Eggs, [theta_e*theta_l*theta_p*F_prod*(1-mu_ad) for k in range(len(Prior_Eggs))])]
                            total_Adults[T] = [i + j + k for i, j, k in zip(total_Adults[T], surviving_Adults[T], new_adults)]
                            total_Eggs[T] = [i + j for i, j in zip(total_Eggs[T], J2)]
                            
                        elif (T_e + T_l) < T < (T_e + T_l + T_p+1):

                             # only the initial adults will be producing offspring at this time.
                            surviving_Adults[T] = [i * j for i, j in zip(total_Adults[T-1], [(1 - k) for k in mortRateAdult1])]
                        
                            proportionPop = convert_to_proportion(surviving_Adults[T])
                            J2 = CRISPR1_1_Equations.CRISPR1_1_Rates(proportionPop, sigma, lamda, delta1, delta2, fitCost, q1, q2,
                                                                     P1, P2, alpha1, alpha2, beta1, beta2, gamma1, gamma2, epsilon1, epsilon2)
                            seeded_Eggs[T] = [i + j for i, j in zip(seeded_Eggs[T], J2)]
                            F_prod = 1
                            for i in range(1,T_l+1):
                                F_prod = F_prod*F(alpha,sum(total_Larvae[T-i]),T_l)
                            new_pupae = [i * j for i, j in zip(seeded_Eggs[T-T_e-T_l], [theta_e*theta_l*F_prod for k in range(len(seeded_Eggs[T-T_e-T_l]))])]
                            surviving_Larvae[T] = [i*j for i, j in zip(total_Larvae[T-1], [(1-mu_l)*F(alpha,sum(total_Larvae[T-1]),T_l) for k in range(len(total_Larvae[T-1]))])]
                            seeded_Larvae[T] = [i * j for i, j in zip(seeded_Eggs[T - T_e], [theta_e for k in range(len(seeded_Eggs[T - T_e]))])]
                            total_Larvae[T] = [i + j + k - q for i, j, k, q in zip(total_Larvae[T], surviving_Larvae[T], seeded_Larvae[T], new_pupae)]
                            # new adults create genotypes of
                            F_prod = 1
                            for i in range(1,T_l+1):
                                F_prod = F_prod*F(alpha,sum(total_Larvae[T-i]),T_l)
                                
                            new_adults = [i * j for i, j in zip(seeded_Eggs[0], [theta_e*theta_l*theta_p*F_prod*(1-mu_ad) for k in range(len(seeded_Eggs[0]))])]
                            total_Adults[T] = [i + j + k for i, j, k in zip(total_Adults[T], surviving_Adults[T], new_adults)]
                            total_Eggs[T] = [i + j for i, j in zip(total_Eggs[T], J2)]
                        
                        else:
                        
                            surviving_Adults[T] = [i * j for i, j in zip(total_Adults[T-1], [(1 - k) for k in mortRateAdult1])]
                            proportionPop = convert_to_proportion(surviving_Adults[T])
                            J2 = CRISPR1_1_Equations.CRISPR1_1_Rates(proportionPop, sigma, lamda, delta1, delta2, fitCost, q1, q2,
                                                                     P1, P2, alpha1, alpha2, beta1, beta2, gamma1, gamma2, epsilon1, epsilon2)
                            seeded_Eggs[T] = [i + j for i, j in zip(seeded_Eggs[T], J2)]
                            # new larvae create genotypes of
                            F_prod = 1
                            for i in range(1,T_l+1):
                                F_prod = F_prod*F(alpha,sum(total_Larvae[T-i]),T_l)
                                
                            new_pupae = [i * j for i, j in zip(seeded_Eggs[T-T_e-T_l], [theta_e*theta_l*F_prod for k in range(len(seeded_Eggs[T-T_e-T_l]))])]
                            surviving_Larvae[T] = [i*j for i, j in zip(total_Larvae[T-1], [(1-mu_l)*F(alpha,sum(total_Larvae[T-1]),T_l) for k in range(len(total_Larvae[T-1]))])]
                            seeded_Larvae[T] = [i * j for i, j in zip(seeded_Eggs[T - T_e], [theta_e for k in range(len(seeded_Eggs[T - T_e]))])]
                            total_Larvae[T] = [i + j + k - q for i, j, k, q in zip(total_Larvae[T], surviving_Larvae[T], seeded_Larvae[T], new_pupae)]
                            # new adults create genotypes of
                            F_prod = 1
                            for i in range(1,T_l+1):
                                F_prod = F_prod*F(alpha,sum(total_Larvae[T-i-T_p]),T_l)
                                
                            new_adults = [i * j for i, j in zip(seeded_Eggs[T-T_e-T_l-T_p], [theta_e*theta_l*theta_p*F_prod*(1-mu_ad) for k in range(len(seeded_Eggs[T-T_e-T_l-T_p]))])]
                            
                            total_Adults[T] = [i + j + k for i, j, k in zip(total_Adults[T], surviving_Adults[T], new_adults)]
                            total_Eggs[T] = [i + j for i, j in zip(total_Eggs[T], J2)]

                        # Store the values of the model into a list
                        total_population_by_genotype[T] = [i + j for i, j in zip(total_population_by_genotype[T], total_Adults[T])]
                        total_population = sum(total_population_by_genotype[T])
                        pop[T] = total_population

                        male_w_count = 2*total_population_by_genotype[T][0] + sum(total_population_by_genotype[T][1:6])
                        male_v_count = total_population_by_genotype[T][1] + 2*total_population_by_genotype[T][6] + sum(total_population_by_genotype[T][7:11])
                        male_u_count = total_population_by_genotype[T][2] + total_population_by_genotype[T][7] + 2*total_population_by_genotype[T][11] + sum(total_population_by_genotype[T][12:15])
                        male_r_count = total_population_by_genotype[T][3] + total_population_by_genotype[T][8] + total_population_by_genotype[T][12] + 2*total_population_by_genotype[T][15] + sum(total_population_by_genotype[T][16:18])
                        male_g_count = total_population_by_genotype[T][4] + total_population_by_genotype[T][9] + total_population_by_genotype[T][13] + total_population_by_genotype[T][16] + 2*total_population_by_genotype[T][18] + total_population_by_genotype[T][19]
                        male_s_count = total_population_by_genotype[T][5] + total_population_by_genotype[T][10] + total_population_by_genotype[T][14] + total_population_by_genotype[T][17] + total_population_by_genotype[T][19] + 2*total_population_by_genotype[T][20]

                        female_w_count = 2*total_population_by_genotype[T][21] + sum(total_population_by_genotype[T][22:27])
                        female_v_count = total_population_by_genotype[T][22] + 2*total_population_by_genotype[T][27] + sum(total_population_by_genotype[T][28:32])
                        female_u_count = total_population_by_genotype[T][23] + total_population_by_genotype[T][28] + 2*total_population_by_genotype[T][32] + sum(total_population_by_genotype[T][33:36])
                        female_r_count = total_population_by_genotype[T][24] + total_population_by_genotype[T][29] + total_population_by_genotype[T][33] + 2*total_population_by_genotype[T][36] + sum(total_population_by_genotype[T][37:39])
                        female_g_count = total_population_by_genotype[T][25] + total_population_by_genotype[T][30] + total_population_by_genotype[T][34] + total_population_by_genotype[T][37] + 2*total_population_by_genotype[T][39] + total_population_by_genotype[T][40]
                        female_s_count = total_population_by_genotype[T][26] + total_population_by_genotype[T][31] + total_population_by_genotype[T][35] + total_population_by_genotype[T][38] + total_population_by_genotype[T][40] + 2*total_population_by_genotype[T][41]

                        w_count = (male_w_count + female_w_count) / (2*total_population)
                        v_count = (male_v_count + female_v_count) / (2*total_population)
                        u_count = (male_u_count + female_u_count) / (2*total_population)
                        r_count = (male_r_count + female_r_count) / (2*total_population)
                        g_count = (male_g_count + female_g_count) / (2*total_population)
                        s_count = (male_s_count + female_s_count) / (2*total_population)

                        output_line = [alpha1, gamma1, alpha2, gamma2, T, w_count, v_count, u_count, r_count, g_count, s_count, total_population]
                        results_temp.append(output_line)

                    print('Run complete in', time.time() - start_time, 'seconds.')
                    
    # START HERE TO OUTPUT THE DATA
    cols = ['male_alpha', 'male_gamma', 'female_alpha', 'female_gamma', 'time', 'w', 'v', 'u', 'r', 'g', 's', 'total_population']
    results_df = pd.DataFrame(results_temp, columns=cols)

    results_df.to_excel('CRISPR1_1 Results.xlsx', index=False)
    
    ## total population plot
    #plt.plot(pop)
    ##plt.yscale("log")
    #plt.ylabel('Total Population', fontsize=8)
    #plt.xlabel('Time Steps', fontsize=8)
    
    ## genotype population plot
    #plots = np.array(total_population_by_genotype)
    #plt.figure()
    #for i in range(2*len(cols)):
    #    plt.plot(plots[:,i])
    #plt.ylabel('Genotype Populations', fontsize=8)
    #plt.xlabel('Time Steps', fontsize=8)

# density dependent process
def F(a,L,T):
    return (a/(a+L))**(1/T)


# converts the list of genotype counts to a proportion by
def convert_to_proportion(K):

    K = [float(i) for i in K]

    halfway_index = len(K)//2
    males = K[:halfway_index]
    females = K[halfway_index:]

    male_total = sum(males)

    if male_total != 0:
        males = [m / male_total for m in males]

    output_list = []
    output_list.extend(males)
    output_list.extend(females)

    return output_list


# generates a list of genotypes based on a list of alleles
def generate_allele_combinations(allele_list):

    # create all the combinations of alleles possible
    allele_1_sequence = []
    allele_2_sequence = []

    allele_dels = allele_list.copy()
    for allele in allele_list:
        for remaining_allele in allele_dels:
            allele_1_sequence.append(allele)
        allele_dels.remove(allele)

    allele_dels = allele_list.copy()
    for allele in allele_list:
        for remaining_allele in allele_dels:
            allele_2_sequence.append(remaining_allele)
        allele_dels.remove(allele)

    # combine the two alleles together into a single genotype
    two_allele_combinations = [i + j for i, j in zip(allele_1_sequence, allele_2_sequence)]

    return two_allele_combinations
