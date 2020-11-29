# historians.py for the historian class to keep track of particle histories

# import reload
from importlib import reload

# import system tools
import os, sys
import json

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

    def __init__(self, electrons=100, status=True, statusii=True):
        """Initialize a historian instance.

        Arguments:
            electrons: int, number of successful electrons
            status: boolean, first slit open?
            statusii: boolean, second slit open?

        Returns:
            None
        """

        # set number of electrons
        self.electrons = electrons

        # set up apparatus
        self.top = 20
        self.bottom = -20
        self.back = -20
        self.divider = 0
        self.screen = 20
        self.source = (-10, 0)

        # configure slits
        self.gap = 1
        self.space = 3
        self.statuses = [status, statusii]
        self.slits = []
        self._configure()

        # set histogram resolution
        self.resolution = 100

        # define probability distribution functions
        self.distribution = lambda x:  2 * cos(pi * x) ** 2
        self.quantile = lambda x: x - 0.5 + sin(2 * pi * x) / (2 * pi)

        # prepare functions for newton-rhapson
        self.zero = lambda x, q: x - q - 0.5 + sin(2 * pi * x) / (2 * pi)
        self.slope = lambda x, q: 1 + cos(2 * pi * x)

        return

    def _bin(self, data, start, finish, resolution):
        """Bin the data from start to finish into a number of bins according to resolution.

        Arguments:
            None

        Returns:
            None
        """

        # create bins
        width = (finish - start) / resolution
        bins = [(start + index * width, start + (index + 1) * width) for index in range(resolution)]

        # get middle points
        middles = [(bin[0] + bin[1]) / 2 for bin in bins]

        # get counts
        counts = [len([datum for datum in data if bin[0] <= datum < bin[1]]) for bin in bins]

        # adjust heights based on area
        area = sum(counts) * width
        heights = [count / area for count in counts]

        return middles, heights

    def _configure(self):
        """Configure the slits based on gap size and space between.

        Arguments:
             None

        Returns:
            None

        Populates:
            self.slits
        """

        # establish first slit location
        top = self.gap + self.space / 2
        bottom = self.space / 2

        # establish second slit location
        topii = -self.space / 2
        bottomii = -self.space / 2 - self.gap

        # set slits
        self.slits = [(top, bottom), (topii, bottomii)]

        return None

    def _crank(self, quantile, guess=1.0, tolerance=1e-14):
        """Crank through Newton Rhapson to solve the quantile equation for a segment length.

        Arguments:
            quantile: float, the quantile to pinpoint
            guess=1.0: float, initial guess
            tolerance=1e-12: tolerance for closesness to zero

        Returns:
            float, better guess for segment length
        """

        # evaluate the function
        evaluation = self.zero(guess, quantile)
        derivative = self.slope(guess, quantile)

        # while the function evaluates to outside the tolerance
        while abs(evaluation) > tolerance:

            # adjust using newton's formula
            guess = guess - evaluation / derivative

            # evaluate the function
            evaluation = self.zero(guess, quantile)
            derivative = self.slope(guess, quantile)

        return guess

    def _cross(self, horizontal, point, pointii):
        """Determine the vertical height at a horizontal between two points.

         Arguments:
             horizontal: float
             point: tuple of floats
             pointii: tuple of floats

        Returns:
            float
        """

        # calculate height
        height = point[1] * (pointii[0] - horizontal) + pointii[1] * (horizontal - point[0]) / (pointii[0] - point[0])

        return height

    def _dump(self, contents, deposit):
        """Dump a dictionary into a json file.

        Arguments:
            contents: dict
            deposit: str, deposit file path

        Returns:
            None
        """

        # dump file
        with open(deposit, 'w') as pointer:

            # dump contents
            json.dump(contents, pointer)

        return None

    def _load(self, path):
        """Load a json file.

        Arguments:
            path: str, file path

        Returns:
            dict
        """

        # open json file
        with open(path, 'r') as pointer:

            # get contents
            contents = json.load(pointer)

        return contents

    def distribute(self):
        """Plot the distribution functions.

        Arguments:
            None

        Returns:
            None
        """

        # get horizontal points
        chunk = 0.01
        horizontals = [number * chunk + 0.5 for number in range(101)]

        # calculate functions
        distributions = numpy.array([self.distribution(horizontal) for horizontal in horizontals])
        quantiles = numpy.array([self.quantile(horizontal) for horizontal in horizontals])

        # set up plot
        pyplot.clf()

        # plot functions
        pyplot.plot(horizontals, distributions, 'b--')
        pyplot.plot(horizontals, quantiles, 'g--')

        # save figure
        pyplot.savefig('distribution.png')

        return None

    def emit(self):
        """Emit electrons from the cathode towards the detector and capture detections

        Arguments:
            None

        Returns:
            None

        Populates:
            self
        """

        # open up tabulated values
        table = self._load('table.json')

        # generate electrons
        while len(self) < self.electrons:

            # begin history at the source
            history = [self.source]

            # random walk until it hits something
            live = True
            keep = False
            while live:

                # pick angle at random
                angle = random() * 2 * pi

                # look up length from the table using a randomly generated quantile
                quantile = str(round(random(), 3))
                length = table[quantile]

                # create new point
                previous = history[-1]
                point = (previous[0] + distance * cos(angle), previous[1] + distance * sin(angle))
                history.append(point)

                # check for hitting top
                if previous[1] <= self.top <= point[1]:

                    # kill it
                    live = False

                # check for hitting bottom
                if point[1] <= self.bottom <= previous[1]:

                    # kill it
                    live = False

                # check for hitting back
                if point[0] <= self.back <= previous[0]:

                    # kill it
                    live = False

                # check for hitting the divider, from either side
                if previous[0] <= self.divider <= point[0] or point[0] <= self.divider <= previous[0]:

                    # kill it
                    live = False

                    # unless it went through the slit
                    cross = self._cross(self.divider, previous, point)
                    if self.slits[0][1] <= cross <= self.slits[0][0] or self.slits[1][1] <= cross <= self.slits[1][0]:

                        # ressurect it
                        live = True

                # check for hitting detector
                if previous[0] <= self.detector <= point[0]:

                    # kill it but keep it
                    live = False
                    keep = True

            # add history to self if keeping
            if keep:

                # append to self
                self.append(history)

                # print status
                if len(self) % 10 == 0:

                    # print status
                    print('{} electrons'.format(len(self), self.electrons))

        # save the histories
        self._dump(self, 'histories.json')

        return None

    def see(self, number=1000):
        """See a number of paths.

        Arguments:
            number: int, max number of histories to plot

        Returns:
            None
        """

        # default number to total length
        number = min([number, len(self)])

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
        for index, history in enumerate(self[-number:]):

            # get coordinates
            xs = numpy.array([point[0] for point in history])
            ys = numpy.array([point[1] for point in history])
            color = colors[index % len(colors)]

            # set opacity
            opacity = 0.01
            for fibonacci in (2, 3, 5, 8, 13):

                # add factor
                opacity += 0.1 * int(index >= number - fibonacci)

            # plot
            pyplot.plot(xs, ys, color, alpha=opacity)

        # plot source
        pyplot.plot(self.source[0], self.source[1], 'yo', markersize=12)
        pyplot.plot(self.source[0], self.source[1], 'go', markersize=9)
        pyplot.plot(self.source[0], self.source[1], 'wo', markersize=4)

        # plot histogram
        chunk = (self.bounds[3] - self.bounds[2]) / self.resolution
        bins = [(self.bounds[2] + index * chunk, self.bounds[2] + (index + 1) * chunk) for index in range(self.resolution)]

        # populate bins
        population = [len([member for member in self if bin[0] < self.cross(self.bounds[1], member[0], member[1]) < bin[1]]) for bin in bins]

        # normalize population
        maximum = max(population)
        population = [entry * 10 / maximum for entry in population]

        # plot it
        for bin, quantity in zip(bins, population):

            # plot the height
            middle = (bin[0] + bin[1]) / 2
            pyplot.plot([self.bounds[1] + 1, self.bounds[1] + 1 + quantity], [middle, middle], 'g-', linewidth=3)
            pyplot.plot([self.bounds[1] + 1, self.bounds[1] + 1 + quantity], [middle, middle], 'y-', linewidth=1)

        # remove ticks
        pyplot.gca().set_xticks([])
        pyplot.gca().set_yticks([])

        # save
        pyplot.savefig('histories.png')

        return

    def tabulate(self, precision=3, tolerance=1e-14):
        """Tabulate values of the quantile integral to be looked up during random walk generation.

        Arguments:
            precision=3: int, precision of table
            tolerance=1e-16: float, tolerance for newton-rhapson convergence

        Returns:
            None
        """

        # generate all quantiles
        quantiles = [float(number) / 10 ** precision for number in range(10 ** precision + 1)]

        # generate all lengths
        lengths = []
        for index, quantile in enumerate(quantiles):

            # print status
            if index % 10 == 0:

                # print status
                print('calculating quantile {} of {}'.format(index, len(quantiles)))

            # calculate length and make a string
            length = self._crank(quantile, guess=1.0, tolerance=tolerance)
            lengths.append(length)

        # create table from rounded quantile strings
        table = {str(round(quantile, precision)): length for quantile, length in zip(quantiles, lengths)}
        self._dump(table, 'table.json')

        return None

    def verify(self, trials=10000):
        """Test the distribution resulting from the inverse quantile function and lookup table.

        Arguments:
            None

        Returns:
            None
        """

        # load in table
        table = self._load('table.json')

        # find all lengths based on random number selection
        lengths = [table[str(round(random(), 3))] for _ in range(trials)]

        # bin them to draw a histogram
        middles, heights = self._bin(lengths, resolution=100, start=0.5, finish=1.5)

        # calculate distribution function
        chunk = 0.01
        horizontals = [number * chunk + 0.5 for number in range(101)]
        distributions = numpy.array([self.distribution(horizontal) for horizontal in horizontals])

        # start plot
        pyplot.clf()

        # plot bars
        for middle, height in zip(middles, heights):

            # plot a bar
            pyplot.plot([middle, middle], [0, height], 'g-', linewidth=1)

        # plot distribution
        pyplot.plot(horizontals, distributions, 'b--')

        # save fig
        pyplot.savefig('verification.png')

        return None

# create instance
historian = Historian(100)
historian.verify()
# historian.emit()
# historian.see()
