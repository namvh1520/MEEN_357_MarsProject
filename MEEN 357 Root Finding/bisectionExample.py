# -*- coding: utf-8 -*-
"""
Created on Thu Sep  8 19:53:43 2022

@author: Christopher Montez
"""

import numpy as np
import matplotlib.pyplot as plt

def basic_bisection(fun, xl, xu, err_max=1e-5, iter_max=1000):
    """
    This function uses bisection to estimate the root of a given function. 

    Parameters
    ----------
    fun : Callable function, dependent on x to be searched
    xl : Scalar, lower bound of root
    xu : Scalar, upper bound of root
    err_max :   Scalar, acceptable approximate relative error.  
                The default is 1e-5.
    iter_max :  Scalar, acceptable number of iterations. The default is 1000.

    Returns
    -------
    root :  Scalar, value of x where function equals zero
    err_est :   Scalar, value of relative error approximation
    numIter :   Scalar, number of iterations used
    exitFlag :  Scalar,
                1 if algorithm terminates properly
                0 if max iteration value is reached
                -1 if invalid brackets given
                -2 if function ever returns with NaN or inf
                
    NOTE :  This function requires that the signs of the endpoints be opposite.

    """
    
    # Checking the inputs are valid
    if not callable(fun):
        raise Exception('The first input must be a callable function.')
        
    if not isinstance(xl, float) and not isinstance(xl, int):
        raise Exception('The second and third inputs must be scalar values.')
        
    if not isinstance(xu, float) and not isinstance(xu, int):
        raise Exception('The second and third inputs must be scalar values.')
        
    if (not isinstance(err_max, float) and not isinstance(err_max, int)) or err_max <=0:
        raise Exception('The fourth and fifth inputs must be positive scalar values.')
        
    if (not isinstance(iter_max, float) and not isinstance(iter_max, int)) or iter_max <=0:
        raise Exception('The fourth and fifth inputs must be positive scalar values.')
        
    # Checking the initial bracket is well-defined
    if xl >= xu:
        raise Exception('The lower bound (second input) must be less than the upper bound (third input).')
    
    # Initializing flag and number of iterations
    isDone = False
    numIter = 0
    
    # Initializing the root estimate and the approximate relative error
    root = np.nan       # Initializing as NaN
    err_est = np.nan
    
    # Checking the endpoints are valid
    fl = fun(xl)
    fu = fun(xu)
    
    if fl * fu > 0:     # Signs are the same, invalid endpoints
        isDone = True
        exitFlag = -1       # Could return here
    elif fl * fu == 0:  # One of the endpoints is already a root
        isDone = True
        err_est = 0.0
        
        if abs(fl) <= abs(fu):  # Finding which one is zero
            root = xl
        else:
            root = xu
        
        # Could return here
        
        
        # No default (i.e., no else) needed
        
    # BEGINNING THE ALGORITHM
    
    while not isDone:
        
        numIter += 1
        
        # Bisecting bracket
        xc = (xu + xl) / 2
        fc = fun(xc)
        
        # Checking the function evaluated at the current guess xc is a valid
        # scalar value. Otherwise, break out of the while loop and return
        # the corresponding exitFlag with NaN as the error estimate and the
        # root estimate. 
        if fc is np.nan:
            isDone = True
            exitFlag = -2
            break
        
        if fl * fc == 0.0:    # We found the root at xc
            isDone = True
            exitFlag = 1
            err_est = 0.0
            root = xc
            break   # Leave the while loop. We could also return here
            
        elif fl * fc > 0:   # fl and fc have same sign, improve lower bound
            xl = xc
        else:               # fu and fc have same sign, improve upper bound
            xu = xc
            
        # NOTE:     We do NOT update the function values fl or fu because 
        #           their sign may change. We just need THEIR INITIAL SIGN to 
        #           figure out if we are updating the lower bound or upper
        #           bound on the root. 
        
        # Updating approximate relative error
        err_est = 100 * abs( (xu - xl) / xc)     # Note the form is different
                                                 # because the width of the
                                                 # bracket tells us how close
                                                 # we are
        
        # NOTE: If you wanted to use the current bracket width and the previous
        #       bracket width to compute an approximate relative error, before
        #       the while loop you would need to initialize a prev_width
        #       (say as np.nan) and a curr_width (initially as xu - xl)
        #       and update these each iteration. Then you can use
        #
        #       err_est = 100 * abs((curr_width - prev_width) / curr_width)
        #
        #       as the approximate relative error for a slightly more
        #       representative measure of error. 
        
        if err_est <= err_max:      # Sufficiently accurate
            isDone = True
            root = xc
            exitFlag = 1
            break       # Could also return instead
            
        if numIter >= iter_max:     # Max number of iterations reached
            isDone = True
            exitFlag = 0
            root = xc   # Current best guess
            break       # Could also return instead
            
    
    # WE HAVE EXITED THE WHILE LOOP
    
    # Notice the indentation level as the same as the keyword while, so this
    # block of text IS NOT IN THE WHILE LOOP
    
    return root, err_est, numIter, exitFlag     # Order matters
    
    
# TESTING THE ALGORITHM

def run():      # Creating a run function to avoid having global variables 
                # (look at the variable explorer after running)

    flagResult = {
        1 : "Terminated normally",
        0 : "Maximum number of iterations reached.",
        -1 : "Invalid bracket given",
        -2 : "Function returned either NaN or Inf"
        }
    
    err_max = 1e-3
    iter_max = 1000
    
    xl = 0
    xu = 2
    
    # Function we wish to find the root(s) of
    def fun(x):
        return x - np.cos(x)
    
    # NOTE: OUTPUT VARIABLES DO NOT NEED TO MATCH THE WORDING IN THE FUNCTION ITSELF
    root, err_est, numIter, exitFlag = basic_bisection(fun, xl, xu, err_max, iter_max)
    
    # FROM WOLFRAM ALPHA, THE ROOT OF x - cos(x) is approximately 0.739085
        
    print("Root: ", root)
    print("Error Estimate: ", err_est)
    print("Number of Iterations: ", numIter)
    print("Exit Flag: ", exitFlag, " -> ", flagResult[exitFlag])
    print("Function Value: ", fun(root))
    
    # VISUALIZATION BY PLOTTING THE FUNCTION
    x = np.linspace(0, 10, 100)
    y = fun(x)
    
    plt.plot(x, y)
    plt.plot(x, 0 * x, color = "black", linestyle = "dashed")
    plt.scatter(root, fun(root), marker="*", color="red", s=100)
    plt.grid()
    plt.rcParams['figure.dpi'] = 600    # Improve resolution of image (try 300)
    plt.show()
    
    
run()
    
    