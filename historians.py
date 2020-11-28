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

    def __init__(self, electrons=100):
        """Initialize a historian instance.

        Arguments:
            None

        Returns:
            None
        """

        # set number of electrons
        self.electrons = electrons

        # set number of deflections per trajectory
        self.deflections = 100

        # initialize detector
        self.bounds = (-20, 20, -20, 20)
        self.source = (-10, 0)
        self.divider = 0
        self.slits = [(-2.5, -1.5), (1.5, 2.5)]

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
        while len(self) < self.electrons:

            # begin history at the source
            history = [self.source]

            # add points to history
            live = True
            keep = False
            while live:

                # pick angle at random
                angle = random() * 2 * pi

                # pick distance as sampling from a normal distribution
                distance = self.sample()

                # create new point
                previous = history[-1]
                point = (previous[0] + distance * cos(angle), previous[1] + distance * sin(angle))
                history.append(point)

                # check for hitting top
                if previous[1] < self.bounds[3] < point[1]:

                    # kill it
                    live = False

                # check for hitting bottom
                if point[1] < self.bounds[2] < previous[1]:

                    # kill it
                    live = False

                # check for hitting back
                if point[0] < self.bounds[0] < previous[0]:

                    # kill it
                    live = False

                # check for hitting detector
                if previous[0] < self.bounds[1] < point[0]:

                    # kill it but keep it
                    live = False
                    keep = True

                # check for hitting the divider
                if previous[0] < self.divider < point[0] or point[0] < self.divider < previous[0]:

                    # kill it
                    live = False

                    # unless it went through the slit
                    cross = previous[1] * (point[0] - self.divider) + point[1] * (self.divider - previous[0]) / (point[0] - previous[0])
                    if self.slits[0][0] < cross < self.slits[0][1] or self.slits[1][0] < cross < self.slits[1][1]:

                        # ressurect it
                        live = True

            # add history to self if keeping
            if keep:

                # append to self
                self.append(history)

                # print status
                if len(self) % 10 == 0:

                    # print status
                    print('{} electrons'.format(len(self), self.electrons))

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

    def see(self, number=None):
        """See a number of paths.

        Arguments:
            number: int, number of histories

        Returns:
            None
        """

        # default number to len
        if not number:

            # default to len
            number = len(self)

        # begin plot
        pyplot.clf()

        # plot detector
        left, right, bottom, top = self.bounds
        pyplot.plot([left, left, right, right, left], [top, bottom, bottom, top, top], 'k-')

        # plot slits
        pyplot.plot([self.divider, self.divider], [bottom, self.slits[0][0]], 'k-')
        pyplot.plot([self.divider, self.divider], [self.slits[0][1], self.slits[1][0]], 'k-')
        pyplot.plot([self.divider, self.divider], [self.slits[1][1], top], 'k-')

        # make colors
        colors = ['r-', 'b-', 'c-', 'm-']

        # plot them
        for index, history in enumerate(self[:number]):

            # get coordinates
            xs = numpy.array([point[0] for point in history])
            ys = numpy.array([point[1] for point in history])
            color = colors[index % len(colors)]

            # set alpha
            alpha = 0.01 + 0.1 * int(index > len(self) - 9) + 0.1 * int(index > len(self) - 6) + 0.1 * int(index > len(self) - 4) + 0.1 * int(index > len(self) - 3) + 0.1 * int(index > len(self) - 2)

            # plot
            pyplot.plot(xs, ys, color, alpha=alpha)

        # plot source
        pyplot.plot(self.source[0], self.source[1], 'yo', markersize=16)
        #pyplot.plot(self.source[0], self.source[1], 'yo', markersize=16)
        pyplot.plot(self.source[0], self.source[1], 'go', markersize=13)
        #pyplot.plot(self.source[0], self.source[1], 'co', markersize=10)
        #pyplot.plot(self.source[0], self.source[1], 'mo', markersize=7)
        pyplot.plot(self.source[0], self.source[1], 'wo', markersize=4)

        # save
        pyplot.savefig('history.png')

        return


# create instance
historian = Historian(1000)
historian.emit()
historian.see()
