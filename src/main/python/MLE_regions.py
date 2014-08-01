#!/usr/bin/python

################################################################################

#Lib,source,work

import sys
import os
import math
from datetime import *

from numpy import *


################################################################################
def matrixrank(A,tol=1e-8):
    s = svd(A,compute_uv=0)
    return sum( where( s>tol, 1, 0 ) )



################################################################################
def fisher( A, T, X, epsilon):
    [m, n] = A.shape
    X_new = X + epsilon
    I = A.T * diag(asarray(1/(A * X_new).T)[0]) * A
    
    
################################################################################
def log_likelihood( A, T, X, P, epsilon):
    X_new = X + epsilon
    common_param_3 = A * X_new;
    f = - sum(common_param_3)+ T * log(common_param_3) + sum(P * X_new)
    
    return [f, common_param_3]


################################################################################
def grad_log_likelihood_i( A, T, X, P, k, common_param_3, epsilon):
    #X_new = X + epsilon
    #Equation: g = - sum(A[:, k]) + sum(  multiply( multiply(asmatrix(T).T,  A[:, k]),  1/(A * X_new) )  )
    g = common_param_1[k] + sum(  divide( common_param_2[k],  (common_param_3) )  ) + P[k][k]


    return g

                             
################################################################################
def fmax_coord_i(A, T, X, P, k, epsilon, step):
    [f, common_param_3] = log_likelihood(A, T, X, P, epsilon)
    g = cmp(grad_log_likelihood_i(A, T, X, P, k, common_param_3, epsilon), 0)
    
    if (g==0):
        return [X ,f, step/2]

    it = 0
    max_it = 10000
    x_old = X[k,0]
    step_old = step

    while (it < max_it):
        X[k] = max(step * g + x_old, 0)
        [new_f, common_param_3] = log_likelihood(A, T, X, P, epsilon)
        if new_f > f:
            return [X , new_f, step]
        if step > epsilon:
            step = step / 2
        else:
            break
        it += 1
        
    if it >= max_it:
        print('warning: maximum iteration exceeded.')

    X[k] = x_old
    return [X, f, step_old]

################################################################################
def fmax_coord( A, T, X, P, n, epsilon):
    [f, common_param_3] = log_likelihood(A, T, X, P, epsilon)
    it = 0
    max_it = 10000
    step = [1e3] * n 

    while (it < max_it):
       
        for i in range(n):
            [X, new_f, new_step] = fmax_coord_i(A, T, X, P, i, epsilon, step[i])
            step[i] = min(new_step * 2, 1e3)
        if abs(new_f - f) < epsilon:
            return [X, new_f]

        f = new_f
        it += 1
        
    if it >= max_it:
        print('warning: maximum iteration exceeded.')
    return [ X, f ]

################################################################################
def solve_likelihood( A, T , P):
    epsilon = 1e-10
    [m, n] = A.shape
    X = asmatrix(ones(n)).T
    [X, f] = fmax_coord(A, T, X, P, n, epsilon)
    #I = fisher(A, T, X, epsilon)
    I = []
    return [ X, f, I]


################################################################################
common_param_1 = []
common_param_2 = []

def compute_common_param( A, N, n):
    global common_param_1
    global common_param_2    
    
    common_param_1 = []
    common_param_2 = []
    for k in range(n):
        common_param_1.append(- sum(A[:, k]) )
        common_param_2.append(multiply(asmatrix(N).T,  A[:, k]))


################################################################################
def iso_opt(A, L, N, NN, P, gene):

    [m, n] = A.shape

    i = 0
    while i < len(L):
        if L[i]<=0:
            #print('warning: L(i)<=0 & N(i)=%d for gene %s.' % (int(N[i]), gene))
            A = delete(A,i,axis=0)
            L = delete(L,i)
            N = delete(N,i)
            continue
        i = i + 1
    
    i = 0
    while i < len(L):
        if all(A[i,]==0):
            #print('warning: A(i,:)==0 & N(i)=%d for gene %s probably new isoform.' % (int(N[i]), gene))
            A = delete(A,i,axis=0)
            L = delete(L,i)
            N = delete(N,i)
            continue
        j = i + 1
        while j < len(L):
            if all(A[i,]==A[j,]):
                A = delete(A,j,axis=0)
                L[i]=L[i]+L[j]
                L = delete(L,j)
                N[i]=N[i]+N[j]
                N = delete(N,j)
            else:
                j = j + 1
        i = i + 1
        
    #if matrixrank(A) < num_isoforms:
        #X_opt = -1
        #print('No MLE computation: num_regions=%d, num_isoforms=%d, num_mapped_reads=%d, rank(A)<n' % (len(L), num_isoforms, int(sum(N))))
        #return [X_opt]
    
    #A = multiply( asmatrix( (NN / 1e6) * (L / 1e3) ).T* asmatrix( ones(n) ), A )
    A = multiply( asmatrix( NN * L * 1e-9 ).T * asmatrix( ones(n) ), A )
    compute_common_param(A, N, n)
    [X_opt, f_max, I]= solve_likelihood(A, N, P)

    return X_opt


################################################################################
def convert_strlist2intlist(strlist):
    
    return [float(item) for item in strlist ]


################################################################################
def parse_penalty_file(penalty_file):
    
    file = open(penalty_file, 'r')
    w = 1
    penalty_dict = {}
    for line in file:
        fields = line.split()
        penalty_dict[int(fields[0])] = float(fields[1]) * w
        
    return penalty_dict


#Main**************************************************************************#
def main():
    # Read input parameters
    input_filename = sys.argv[1]
    output_filename = sys.argv[2]
    if (len(sys.argv) > 3):
        penalty_file = sys.argv[3]
        penalty_dict = parse_penalty_file(penalty_file)
        add_penalty = True
    else:
        add_penalty = False
        
        
    input_file = open(input_filename,"r")
    output_file = open(output_filename,"w")
    
    t0=datetime.now()
    header = input_file.readline()
    NN = float(header.split()[0])

    while True:
        
        gene_line = input_file.readline().strip()
        if gene_line == "": 
            break
        gene_line_list = gene_line.split()
        gname = gene_line_list[0]
        chrnum = gene_line_list[2]
    
        num_isoforms = int (gene_line_list[1])
        isoform_names = input_file.readline().strip().split()
        isoform_known = [int(i) for i in input_file.readline().strip().split()] 
        isoform_lengths = [int(l) for l in input_file.readline().strip().split()]
    
        input_file.readline()
        input_file.readline()   
        input_file.readline()
        
        A = []
        for i in range(num_isoforms):
            A.append( convert_strlist2intlist ( input_file.readline().strip().split() ) )
        A = asmatrix( array(A) ).T
        
        L = array( convert_strlist2intlist ( input_file.readline().strip().split() ) )
        N = array( convert_strlist2intlist ( input_file.readline().strip().split() ) )
        len_ljust = max((len(gname)/10 + 1) * 10, 20)
        output_file.write(gname.ljust(len_ljust) + chrnum.ljust(10) + \
                          str(len(L)).ljust(10) + str(num_isoforms).ljust(10) + \
                          str(int(sum(N))).ljust(10) + '\n')
        
        P = [0] * num_isoforms
        if (add_penalty):
            for i in range(num_isoforms):
                if (isoform_known[i] == 0):
                    P[i] = penalty_dict[isoform_lengths[i]]
        P = diag(array(P), 0)
        
        X_opt = -1
        if (num_isoforms > 1) and (len(L) > 1) and (sum(N) >= 10):
            X_opt = iso_opt(A, L, N, NN, P, gname)
        #else:
            #print('No MLE computation: num_regions=%d, num_isoforms=%d, num_mapped_reads=%d' % (len(L), num_isoforms, int(sum(N))))
        
        str_print = ""
        for i in range(num_isoforms):
            len_ljust = max((len(isoform_names[i])/10 + 1) * 10, 20)
            str_print += (isoform_names[i].ljust(len_ljust))
        output_file.write(str_print + ('Total').ljust(20) + '\n')
    
        if any(X_opt < 0):    
            N = sum(N)
            L = sum(L)
            #[X_opt] = [N / ((L) * 1e-3) / (NN * 1e-6)]
            X_opt = N / (L) / NN * 1e9
            if (num_isoforms == 1):
                X_opt_list = [[X_opt]]
            else:
                X_opt_list = ['NA'] * num_isoforms
            total_exp_level = X_opt
            #print 'Total Expression: ' + str(X_opt)
        else:
            #print('MLE computation: num_regions=%d, num_isoforms=%d, num_mapped_reads=%d' % (len(L), num_isoforms, int(sum(N))))
            X_opt_list = X_opt.tolist() 
            total_exp_level = sum(X_opt_list)
        str_print = ""

        for i in range(num_isoforms):
            len_ljust = max((len(isoform_names[i])/10 + 1) * 10, 20)
            if (X_opt_list[i] == 'NA'):
                str_print += (str(X_opt_list[i]).ljust(len_ljust))
            else:
                str_print += (str(round(X_opt_list[i][0], 4)).ljust(len_ljust))
        output_file.write(str_print + str(round(total_exp_level, 4)).ljust(20))
        output_file.write('\n')
    
    input_file.close()
    output_file.close()
    
    t1=datetime.now()
    print "get RNA info running time: " + str(t1-t0)



if __name__ == '__main__':
    main()