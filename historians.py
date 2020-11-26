# historians.py for the historian class to keep track of particle histories

# import reload
from importlib import reload

# import system tools
import os, sys

# import random
from random import random

# import trig functions
from math import sin, cos, pi

# import numpy
import numpy

# import matplotlib for plots
from matplotlib import pyplot
from matplotlib import gridspec
from matplotlib import colorbar
from matplotlib import colors as Colors
from matplotlib import style as Style
Style.use('fast')


# class historian
class Historian(list):
    """Class Historian to keep track of randomly generated electron histories.

     Inherits from:
        list
    """

    def __init__(self):
        """Initialize a historian instance.

        Arguments:
            None

        Returns:
            None
        """

        # set number of electrons
        self.electrons = 1000

        # set number of deflections per trajectory
        self.deflections = 100

        return

    def emit(self):
        """Emit electrons from the cathode towards to detector.

        Arguments:
            None

        Returns:
            None

        Populates:
            self
        """

        # generate electrons
        for trial in range(self.electrons):

            # print status
            if trial % 100 == 0:

                # print status
                print('trial {} of {}...'.format(trial, self.electrons))

            # begin history at the origin
            history = [(0.0, 0.0)]

            # add points to history
            for _ in range(self.deflections):

                # pick angle at random
                angle = random() * 2 * pi

                # pick distance as sampling from a normal distribution
                distance = self.sample()

                # create new point
                previous = history[-1]
                point = (previous[0] + distance * cos(angle), previous[1] + distance * sin(angle))
                history.append(point)

            # add history to self
            self.append(history)

        return None

    def sample(self):
        """Sample a distance from the sin^2 distribution, using the newton rhapson method.

        Arguments:
            None

        Returns:
            float
        """

        # consider a normalesque probability density function p = (sin (pi x / 2))^2
        # it has it's maximum at 1 and tapers off to 0 at 0 and 2
        # the integral of this function describes a cumulative distribution c = (pi x - sin(pi x)) / 2 pi

        # consider a number u, sampled at random from a uniform distribution from 0 to 1.  Setting this equal to the
        # cumulative distribution function and solving for x will be the distance leading to the sin^2 distribution

        # this is not solvable exactly, however, so newton's method will be used:
        # x_n+1 = x_n - f(x_n) / f'(x_n), with an initial guess of 1
        # u = (pi x - sin(pi x)) / 2 pi
        # f(x) = pi x - sin(pi x)) - 2 pi u = 0
        # f'(x) = pi - pi cos(pi x))

        # get a random number from 0 to 1
        number = random()

        # define the function to zero, and its derivative
        zero = lambda x: pi * x - sin(pi * x) - 2 * pi * number
        slope = lambda x: pi - pi * cos(pi * x)

        # define the tolerance
        tolerance = 1e-12

        # begin with initial guess 1
        guess = 1.0

        # evaluate the function
        evaluation = zero(guess)
        derivative = slope(guess)

        # while the function evaluates to outside the tolerance
        while abs(evaluation) > tolerance:

            # adjust using newton's formula
            guess = guess - evaluation / derivative

            # evaluate the function
            evaluation = zero(guess)
            derivative = slope(guess)

        return guess

    def see(self, number):
        """See a number of paths.

        Arguments:
            number: int, number of histories

        Returns:
            None
        """

        # begin plot
        pyplot.clf()

        # make colors
        colors = ['k--', 'r--', 'g--', 'b--', 'c--', 'm--']

        # plot them
        for index, history in enumerate(self[:number]):

            # get coordinates
            xs = numpy.array([point[0] for point in history])
            ys = numpy.array([point[1] for point in history])
            color = colors[index % len(colors)]

            # plot
            pyplot.plot(xs, ys, color)

        # save
        pyplot.savefig('history.png')

        return


# create instance
historian = Historian()
historian.emit()
historian.see(1000)
