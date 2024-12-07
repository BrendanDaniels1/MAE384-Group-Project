SIR model Project Overview:

Parts 1 & 2 : Ian
Problems one and two required using the SIR model to analyze disease dynamics for measles, influenza, and COVID-19. The three differential equations represented the populations of susceptible (not yet infected), infected (actively infected), and recovered individuals (those reintroduced into the population, assuming no deaths).
Using the Runge-Kutta method with daily step sizes, we approximated how quickly individuals became infected and recovered for each disease. For part of the problem, we used every other data point from these results and interpolated the missing values using a 4th-degree polynomial. We then compared the interpolated values to the original data and calculated the percent error, achieving very accurate results.

Part 3: Isela Gonzales
This part involves using the least squares regression method in order to estimate the model parameter, k, as well as the initial condition, I0,
of the infected population. These estimates, with different data set sizes, were then compared to the true data from part 1 
Discussion:
Based on the results, the data from merely 10 days of running the model
deemed to be more accurate in comparison to 30 days of data. This could
be due to the values needed being at the beginning stages of the model
and the linearization of this model makes the more accurate numbers in
the beginning.

Part 4: Brendon
