# -*- coding: utf-8 -*-
"""
Created on Tue Sep 13 11:39:10 2022

@author: Nam Vo
"""
import numpy as np
import matplotlib.pyplot as plt


def falsepos(fun, lb, ub, err_max = 1e-5, iter_max =1000):
    
    #Checking if input are valid
    if not callable(fun):
        raise Exception('The first input must be a callable function.')
        
    if not isinstance(lb, float) and not isinstance(lb, int):
        raise Exception('The second and third inputs must be scalar values.')
        
    if not isinstance(ub, float) and not isinstance(ub, int):
        raise Exception('The second and third inputs must be scalar values.')
        
    if (not isinstance(err_max, float) and not isinstance(err_max, int)) or err_max <=0:
        raise Exception('The fourth and fifth inputs must be positive scalar values.')
        
    if (not isinstance(iter_max, float) and not isinstance(iter_max, int)) or iter_max <=0:
        raise Exception('The fourth and fifth inputs must be positive scalar values.')
    
    if lb >= ub:
        raise Exception('The lower bound (second input) must be less than the upperbound(third input)')

    isDone = False
    numIter = 0
    
    root = np.nan
    err_est = np.nan
    
    fl = fun(lb)
    fu = fun(ub)
    
    if fl * fu > 0:
        isDone = True
        exitFlag = -1
    elif fl * fu == 0:
        isDone = True
        err_est = 0.0
        exitFlag = 1
        if abs(fl) <= abs(fu):
            root = lb
        else:
            root = ub
            
    while not isDone:
        
        numIter += 1
        
        #Getting function value
        fl = fun(lb)
        fu = fun(ub)
        
        #Finding xr
        xr = ub - (fu * (lb - ub) /(fl - fu))
        
        fc = fun(xr)
        
        #Checking if return is valid
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
        elif fl * fc > 0:
            lb = xr
        else: 
            ub = xr
            
        err_est = 100 * abs((ub - lb)/xr)
            
        if err_est <= err_max:
            isDone = True
            root = xr
            exitFlag = 1
            
        if numIter >= iter_max:
            isDone = True
            exitFlag = 0
            root = xr
            break
        
    return root, err_est, numIter, exitFlag


def secant(fun, ini_guess, err_max = 1e-5, iter_max =1000):
    
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
            
        err_est = 100 * abs((ini_guess - guess2)/ini_guess)
            
        if err_est <= err_max:
            isDone = True
            root = xr
            exitFlag = 1
            
        if numIter >= iter_max:
            isDone = True
            exitFlag = 0
            root = xr
            break
        
    return root, err_est, numIter, exitFlag


def run():     

    flagResult = {
        1 : "Terminated normally",
        0 : "Maximum number of iterations reached.",
        -1 : "Invalid bracket given",
        -2 : "Function returned either NaN or Inf"
        }
    
    err_max = 10e-100
    iter_max = 1000
    
    xl = 0
    xu = 30

    def funfun(x):
        return (x**3) / 88
    
    def fun1(x):
        return x * np.sin(x) + 3 *np.cos(x) - x
    
    def fun2(x):
        return x * (np.sin(x) - x * np.cos(x))
    
    def fun3(x):
        return (x**3 - 2 * x**2 + 5 * x -25)/40
    
    root, err_est, numIter, exitFlag = secant(funfun, xl,err_max, iter_max)
    
    # Check with Wolfram Alpha or other applications for correct root
        
    print("Root: ", root)
    print("Error Estimate: ", err_est)
    print("Number of Iterations: ", numIter)
    print("Exit Flag: ", exitFlag, " -> ", flagResult[exitFlag])
    print("Function Value: ", funfun(root))
    
    x = np.linspace(-6, 30, 100)
    y = funfun(x)
    
    plt.plot(x, y)
    plt.plot(x, 0 * x, color = "black", linestyle = "dashed")
    plt.scatter(root, funfun(root), marker="*", color="red", s=100)
    plt.grid()
    plt.rcParams['figure.dpi'] = 300
    plt.show()
    
run()