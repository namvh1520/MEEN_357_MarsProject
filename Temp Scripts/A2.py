
import math
import numpy as np
def falsepos(fun, xl, xu, err_max = 1e-5, iter_max = 1000):
    
    #check for exceptions
    if not callable(fun):
        raise Exception('The first input must be a function.')
        
    if not isinstance(xl, float) and not isinstance(xl, int):
        raise Exception('The second and third inputs must be scalar values.')
        
    if not isinstance(xu, float) and not isinstance(xu, int):
        raise Exception('The second and third inputs must be scalar values.')
        
    if (not isinstance(err_max, float) and not isinstance(err_max, int)) or err_max<=0:
        raise Exception('The fourth and fifth inputs must be positive scalar values.')
        
    if (not isinstance(iter_max, float) and not isinstance(iter_max, int)) or iter_max <=0:
        raise Exception('The fourth and fifth inputs must be positive scalar values.')
        

    if xl >= xu:
        raise Exception('The lower bound (second input) must be less than the upperbound (third input).')    
    
    #first assignment of variables
    root = math.nan
    err = math.nan
    numIter = 0
    exitFlag = math.nan
    
    
    isDone = False
    
  #if the bounds contain no root
    '''
    if fun(xl) * fun(xu) > 0:
        isDone = True
        exitFlag = -1
        raise Exception('There is no root in stated domain')
  
    #if one of the bounds is the root
    if fun(xl) * fun(xu) == 0:
        isDone = True
        err = 0
        if abs(fun(xl)) < abs(fun(xu)):
            root = xl
        else:
            root = xu
    '''
    #algorithm
    while not isDone:
        
        numIter += 1
        
        root_x = xu - (fun(xu)*(xl-xu))/(fun(xl)-fun(xu))
        
        err = (abs(root_x - xl) / abs(root_x))
        
        if fun(root_x) is math.nan:
            isDone = True
            exitFlag = -2
            break
        
        if fun(root_x) == 0:
            exitFlag = 1
            err = 0
            isDone = True
            root = root_x
            break
            
        # the guess takes the place of the lower bound
        elif root_x * xl >= 0:
            xl = root_x
        # the guess takes the place of the upper bound
        else:
            xu = root_x
    
        
        if err <= err_max:
            exitFlag = 1
            isDone = True
            root = root_x
            break
        
        if numIter >= iter_max:
            exitFlag = 0
            isDone = True
            root = root_x
            print('No root was found')
            break
    
    
    
    
    return root, err, numIter, exitFlag



def secant(fun, ini_guess, err_max = 10e-100, iter_max =10000):
    
    if not callable(fun):
        raise Exception('The first input must be a callable function.')
        
    if not isinstance(ini_guess, float) and not isinstance(ini_guess, int):
        raise Exception('The second inputs must be scalar values.')
                
    if (not isinstance(err_max, float) and not isinstance(err_max, int)) or err_max <=0:
        raise Exception('The fourth and fifth inputs must be positive scalar values.')
        
    if (not isinstance(iter_max, float) and not isinstance(iter_max, int)) or iter_max <=0:
        raise Exception('The fourth and fifth inputs must be positive scalar values.')
    
    guess2 = ini_guess + .001
    isDone = False
    numIter = 0
    
    root = np.nan
    err_est = np.nan
    
    fl = fun(ini_guess)
    fu = fun(guess2)
    
    if fl * fu == 0:
        isDone = True
        err_est = 0.0
        exitFlag = 1
        if abs(fl) <= abs(fu):
            root = ini_guess
        else:
            root = guess2
        
    while not isDone: 
        numIter += 1
        
        fl = fun(ini_guess)
        fu = fun(guess2)
        
        xr = ini_guess - (fun(ini_guess) * (guess2 - ini_guess)/(fun(guess2) - fun(ini_guess)))
        
        fc = fun(xr)
        
        if fc is np.nan:
            isDone = True
            exitFlag = -2
            break
        
        
        if fl * fc == 0.0:
            isDone = True
            exitFlag = 1
            err_est = 0.0
            root = xr
            break
        else: 
            guess2 = ini_guess
            ini_guess = xr
            
        err_est = 100 * abs((xr - ini_guess)/xr)
            
        if err_est <= err_max:
            isDone = True
            root = xr
            exitFlag = 1
            
        if numIter >= iter_max:
            isDone = True
            exitFlag = 0
            root = xr
            break
        
    return root#, err_est, numIter, exitFlag

    
    
    
    
#code for testing either method

