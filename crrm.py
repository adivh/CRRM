##################################################
### Please use python3.x to avoid inaccuracies ###
### caused by integer division.   Thank you!   ###
##################################################

from scipy.stats import norm
from math import exp
from math import log
from math import nan
from math import sqrt
import matplotlib.pyplot as plt
import sys

if sys.version_info[0] < 3:
    print("Please use Python 3.x")
    print("exiting..")
    sys.exit()

# parameters
# V       verbose output
# PLOT    plot function n: [1..PLOT] -> price: []
#         set to 0 for simple depth-based calculation
#         set to value > 0 for a plot with n in [1,value]
# PCS.    0 = Put, 1 = Call, 2 = Sprint
# T...    expiration time in years
# S...    stock price
# E...    strike price
# BAS.    used for sprint: yield above BAS is increased by factor mult up until CAP
# CAP.    used for sprint: derivative does not yield additional returns when S > CAP
# MUL.    multiplier for sprint: derivative yields MUL times more returns for range E < S < CAP
#         sprint price is = {S <= E: S, E < S <= CAP: E + MUL * (S - E), CAP < S: E + 2 * (CAP - E)}
# q...    dividend yield NOT YET IMPLEMENTED
# n...    height of the binomial tree
# r...    interest rate
# sigma.. volatility
# mu..    expected return of log norm

V = 1
PLOT = 000
PCS = 0
T = 4/12            # 3 months remaining
S = 119.1135
E = 120.0000
BAS = 126.0000
CAP = 130.0000
MUL = 2
q = 0               # not implemented
n = 100000
r = -.04            # 3 months EURIBOR
sigma = .3086
mu = -0.0

# catch invalid PCS value
# PCS must be {1,2,3}

if PCS < 0 or PCS > 2:
    print("PCS={} is not defined. Use PCS = 1,2,3".format(PCS))
    sys.exit()

# calculate option value at expiration

def intrinsic_derivative_value(s):
    if PCS == 0:
        # put option
        return max(E - s, 0)
    
    elif PCS == 1:
        # call option
        return max(s - E, 0)

    elif PCS == 2:
        # sprint certificate
        if s <= BAS:
            return s

        elif s <= CAP:
            return E + MUL * (s - BAS)

        else:
            return CAP * 2 - BAS
    else:
        # PCS not defined
        print('this should never be reached..')
        return nan
        

# used for recursive calculation in the binomial model

def node_value(s1, s2, p, exprdt):
    pr = 1 - p

    # return discounted expected return
    
    return (p * s1 + pr * s2) * exprdt
    
# Calculate Cox-Ross-Rubinstein Binomial Model
# Return current price based on model

def crrm(n):

    deltaT = T / n

    # u... expected change for state up
    # d... expected change for state down
    # p... probability for up; 1-p probability for down
    #
    #      expected change is change up   + change down     - original value
    #                         [p * S * u] + [(1-p) * S * d] - [S]

    u = exp(sigma * sqrt(deltaT))
    d = 1 / u

    # CRR probability calculation

    p = .5 + mu * sqrt(deltaT) / 2 / sigma

    # precomputation of discount factor

    exprdt = exp(-r * deltaT)

    # catch errors
    # this should only occure for wierd values in mu or sigma

    if p <= 0 or p >= 1:
        print('probability for up-state not in (0,1)..')
        print('mu={}, sigma={}'.format(mu, sigma))
        print('exiting..')
        sys.exit()

    # print values for fun

    if V > 0:
        print('u={}, d={}, p={}'.format(u, d, p))

    # initialize 2d list in imperative language style because I don't know any better

    prices = [[nan] * (n + 1), [nan] * n]

    # Algorithm:
    # prune the tree as soon as the leaves become obsolete
    # this reduces memory complexit from O(n^2) to O(n)
    #
    # iteration 1                                      iteration 2                                       --->    iteration n                                   
    # calculate parent         write back              calculate parent          write back              --->    calculate parent          write back          
    #
    # t=n          t=n-1       t=n         t=n-1       t=n-1         t=n-2       t=n-1       t=n-2       --->    t=n-1         t=n-2       t=n-1       t=n-2   
    #                                                                                                    --->                             !!!!!!!!             
    # [0][0] --+-> [1][0]      [0][0] <--- [1][0]      [0][0] --+--> [1][0]      [0][0] <--- [1][0]      --->    [0][0] --+--> [1][0]      RESULT <--- [1][0]  
    #          |                                                |                                        --->             |               !!!!!!!!             
    # [0][1] --/   [1][1]      [0][1]      [1][1]      [0][1] --/    [1][1]      [0][1]      [1][1]      --->    [0][1] --/    [1][1]      [0][1]      [1][1]  
    # [0][3]       [1][2]      [0][3]      [1][2]      [0][3]        [1][2]      [0][3]      [1][2]      --->    [0][3]        [1][2]      [0][3]      [1][2]  
    # [0][4]       [1][3]      [0][4]      [1][3]      [0][4]        [1][3]      [0][4]      [1][3]      --->    [0][4]        [1][3]      [0][4]      [1][3]  
    #   ..           ..          ..          ..          ..            ..          ..          ..        --->      ..            ..          ..          ..    
    # [0][n-2]     [1][n-2]    [0][n-2]    [1][n-2]    [0][n-2]  +-> [1][n-2]    [0][n-2] <- [1][n-2]    --->    [0][n-2]      [1][n-2]    [0][n-2]    [1][n-2]
    #                                                            |                                       --->                                                  
    # [0][n-1] +-> [1][n-1]    [0][n-1] <- [1][n-1]    [0][n-1] -/   [1][n-1]    [0][n-1]    [1][n-1]    --->    [0][n-1]      [1][n-1]    [0][n-1]    [1][n-1]
    #          |                                                                                         --->                                                  
    # [0][n] --/               [0][n]                  [0][n]                    [0][n]                  --->    [0][n]                    [0][n]              

    # calculate option prices for each state at expiration
    # save result in leaves

    for i in range(0, n + 1):
        prices[0][i] = intrinsic_derivative_value(S * u ** (n - 2 * i))

    # traverse tree from leaves to root
    # calculate node values based on the child node values

    for i in range(1, n):

        # output progress

        if V > 0 and i % 250 == 0:
            print('{:5d}'.format(i))

        # calculate node value based on both child node values and save in current node
        # then write back result in what used to represent the upper child node

        for j in range(0, n - i):
            prices[1][j] = node_value(prices[0][j], prices[0][j + 1], p, exprdt)
            prices[0][j] = prices[1][j]

    # the root is the current price
    # the root is now stored in [0][0]

    return prices[0][0]

# Caluclate Black-Scholes Model
# Return current price based on model

def bs():
    d1 = log(S/E) + (r + sigma * sigma / 2) * T / sigma / sqrt(T)
    d2 = d1 - sigma * sqrt(T)

    if PCS == 1:
        # call option
        return S * norm.cdf(d1) - E * exp(-r * T) * norm.cdf(d2)

    elif PCS == 0:
        # put option
        return E * exp(-r * T) * norm.cdf(-d2) - S * norm.cdf(-d1)

    else:
        # this should never be reached
        print('PCS={}, PCS must be 0 or 1 to use Black-Scholes')
        return nan

# no plot
# if option -> calculate CRRM and BS
# if sprint -> calculate CRRM

if PLOT <= 0:
    crrm = crrm(n)
    print('Aktueller Preis nach CRRM mit n={}: {}'.format(n, crrm))

    if PCS == 0 or PCS == 1:
        bs = bs()
        print('Aktueller Preis nach Black-Scholes: {}'.format(bs))

# yes plot
# if option -> calculate CRRM and BS
# if sprint -> calculate CRRM

else:
    # collect price data for every depth from n to PLOT

    price_per_depth = []

    for i in range(1,PLOT + 1):

        if i % 100 == 0:
            print(i)

        price_per_depth.append(crrm(i))

    # x axis 1 to PLOT
    # will be mapped to price per depth

    x = list(range(1,PLOT + 1))

    # plot price calculated with CRRM

    plt.plot(x, price_per_depth, label='CRRM price (mu = 0)')

    # plot price caluclated with BS and actual price for comparison

    if PCS == 0 or PCS == 1:
        plt.plot(x, [bs()] * PLOT, label='Black-Scholes price')
        plt.plot(x, [9.1750] * PLOT, label='Actual price')
    
    # apply logarithmical scale to depth values on the x axis

    plt.xscale("log")

    # make the plot prettier

    plt.xlabel('depth of binomial tree')
    plt.ylabel('price [EUR]')
    plt.legend()

    # render the plot

    plt.show()
