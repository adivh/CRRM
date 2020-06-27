##################################################
### Please use python3.x to avoid inaccuracies ###
### caused by integer division.   Thank you!   ###
##################################################

from scipy.stats import norm
from math import exp
from math import log
from math import sqrt
import matplotlib.pyplot as plt
import sys

if sys.version_info[0] < 3:
    print("Please use Python 3.x")
    print("exiting..")
    sys.exit()

# parameters
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
# q...    dividend yield
# n...    height of the binomial tree
# r...    interest rate
# sigma.. volatility
# mu..    expected return of log norm

PLOT = 000
PCS = 2
T = 4/12            # 3 months remaining
S = 119.1135
E = 120.0000
BAS = 126.0000
CAP = 130.0000
MUL = 2
q = 0
n = 1000
r = -.04            # 3 months EURIBOR
sigma = .3086
mu = -0.0

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
        # CP not defined

        print("PCS={} is not defined. Use PCS = 1,2,3".format(PCS))
        sys.exit()

# used for recursive calculation in the binomial model
def node_value(s1, s2, p, r, dt):
    pr = 1 - p

    # return discounted expected return

    return (p * s1 + pr * s2) * exp(-r * dt)
    

# calculate option prices for each state at expiration

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

    # catch errors
    # this should only occure for wierd values in mu or sigma

    if p <= 0 or p >= 1:
        print('probability for up-state not in (0,1)..')

    # print values for fun

    if PLOT <= 0:
        print('u={}, d={}, p={}'.format(u, d, p))

    # initialize 2d list in imperative language style because I don't know any better

    prices = [[]]

    # STRUCTURE:
    # t=n       t=n-1   ..  t=1         t=1
    # [0][0]    [1][0]      [n-1][0]    [n][0]
    # [0][1]    [1][1]      [n-1][1]
    # [0][3]    [1][2]
    # [0][4]    [1][3]
    # [0][5]    [1][4]
    # [0][6]    [1][5]
    # [0][7]

    for i in range(0, n + 1):
        prices[0].append(intrinsic_derivative_value(S * u ** (n - 2 * i)))

    for i in range(1, n):
        prices.append([])
        for j in range(0, n - i):
            prices[i].append(node_value(prices[i - 1][j], prices[i - 1][j + 1], p, r, deltaT))

    if PLOT <= 0:
        print('Aktueller Preis nach CRRM mit n={}: {}'.format(n, prices[-1][0]))
    else:
        return prices[-1][0]

# Black-Scholes

def bs():
    d1 = log(S/E) + (r + sigma * sigma / 2) * T / sigma / sqrt(T)
    d2 = d1 - sigma * sqrt(T)

    if PCS > 0:
        bs_price = S * norm.cdf(d1) - E * exp(-r * T) * norm.cdf(d2)

    else:
        bs_price = E * exp(-r * T) * norm.cdf(-d2) - S * norm.cdf(-d1)

    print('Aktueller Preis nach Black-Scholes: {}'.format(bs_price))

# Run price calculation

if PLOT <= 0:
    crrm(n)
    if PCS == 0 or PCS == 1:
        bs()

else:
    price_per_depth = []
    for i in range(1,PLOT + 1):
        price_per_depth.append(crrm(i))

    plt.plot(list(range(1,PLOT + 1)), price_per_depth)
    plt.xscale("log")
    plt.xlabel('depth of binomial tree')
    plt.ylabel('price [EUR]')
    plt.show()
