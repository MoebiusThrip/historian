# historians.py for the historian class to keep track of particle histories

# import reload
from importlib import reload

# import system tools
import os, sys
import json

# import time
from time import time

# import random

# import trig functions
from math import sin, cos, pi, exp, log

# import numpy
import numpy
from numpy.random import random
from numpy import sin, cos, pi, exp, log

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

        # current time
        self.time = time()

        # set number of electrons
        self.electrons = electrons

        # set up apparatus
        self.top = 30
        self.bottom = -30
        self.back = -20
        self.divider = 0
        self.screen = 30
        self.source = (-10, 0)

        # configure slits
        self.gap = 1
        self.space = 7
        self.statuses = [status, statusii]
        self.slits = []
        self._configure()

        # set histogram resolution
        self.resolution = 10000

        # set tabulation precision
        self.precision = 3

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

    def _build(self):
        """Plot the detector apparatus in the plot.

        Arguments:
            None

        Returns:
            None
        """

        # plot detector
        left = self.back
        right = self.screen
        top = self.top
        bottom = self.bottom
        divider = self.divider
        pyplot.plot([left, left, right, right, left], [top, bottom, bottom, top, top], 'k-', linewidth=3)
        pyplot.plot([left, left, right, right, left], [top, bottom, bottom, top, top], 'g-', linewidth=1)

        # plot slits
        pyplot.plot([divider, divider], [bottom, self.slits[1][1]], 'k-', linewidth=3)
        pyplot.plot([divider, divider], [self.slits[0][1], self.slits[1][0]], 'k-', linewidth=3)
        pyplot.plot([divider, divider], [self.slits[0][0], top], 'k-', linewidth=3)

        # plot slits
        pyplot.plot([divider, divider], [bottom, self.slits[1][1]], 'g-', linewidth=1)
        pyplot.plot([divider, divider], [self.slits[0][1], self.slits[1][0]], 'g-', linewidth=1)
        pyplot.plot([divider, divider], [self.slits[0][0], top], 'g-', linewidth=1)

        # add light to screen
        pyplot.plot([right, right], [top, bottom], 'g-', linewidth=2)

        return None

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
        print('dumping into {}...'.format(deposit))
        with open(deposit, 'w') as pointer:

            # dump contents
            json.dump(contents, pointer)

        return None

    def _graph(self):
        """Plot the histograph.

        Arguments:
            None

        Returns:
            None
        """

        # status
        print('plotting histogram...')

        # plot histogram
        chunk = (self.top - self.bottom) / self.resolution
        bins = [(self.bottom + index * chunk, self.bottom + (index + 1) * chunk) for index in range(self.resolution)]

        # populate bins
        population = [len([member for member in self if bin[0] < self._cross(self.screen, member[0], member[1]) < bin[1]]) for bin in bins]

        # normalize population
        maximum = max(population)
        population = [entry * 10 / maximum for entry in population]

        # plot it
        for bin, quantity in zip(bins, population):

            # plot the height
            middle = (bin[0] + bin[1]) / 2
            pyplot.plot([self.screen + 1, self.screen + 1 + quantity], [middle, middle], 'g-', linewidth=0.5)

        return None

    def _ignite(self):
        """Plot the source.

        Arguments:
            None

        Returns:
            None
        """

        # plot source
        pyplot.plot(self.source[0], self.source[1], 'yo', markersize=12)
        pyplot.plot(self.source[0], self.source[1], 'go', markersize=9)
        pyplot.plot(self.source[0], self.source[1], 'wo', markersize=4)

        return None

    def _load(self, path):
        """Load a json file.

        Arguments:
            path: str, file path

        Returns:
            dict
        """

        # open json file
        print('loading {}...'.format(path))
        with open(path, 'r') as pointer:

            # get contents
            contents = json.load(pointer)

        return contents

    def _time(self, message):
        """Start timing a block of code, and print results with a message.

        Arguments:
            message: str

        Returns:
            None
        """

        # get final time
        final = time()

        # calculate duration and reset time
        duration = final - self.time
        self.time = final

        # print duration
        print('took {} seconds.'.format(duration))

        # begin new block
        print(message)

        return None

    def _trace(self, trajectories):
        """Trace trajectories on the plot.

        Arguments:
            trajectories: tuples

        Returns:
            None
        """

        # status
        print('plotting trajectories...')

        # plot the trajectories
        for history, color, opacity in trajectories:

            # get horizontals and verticals
            horizontals = numpy.array([point[0] for point in history])
            verticals = numpy.array([point[1] for point in history])

            # plot
            pyplot.plot(horizontals, verticals, color, alpha=opacity, linewidth=0.5)

        return None

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

        # generate set of random floats
        randoms = random(size=1000000).tolist()

        # generate electrons
        count = 0
        start = time()
        while len(self) < self.electrons:

            # begin history at the source
            history = [self.source]
            count += 1

            self._time('starting...')

            # random walk until it hits something
            live = True
            keep = False
            while live:

                self._time('starting electron...')

                # generate more random numbers
                if len(randoms) < 10:

                    # generate more random numbers
                    randoms = random(size=1000000).tolist() + randoms

                # pick angle at random
                angle = randoms.pop() * 2 * pi

                # look up length from the table using a randomly generated quantile
                index = int(randoms.pop() * 10 ** self.precision)
                length = table[index]

                self._time('made randoms')

                # create new point
                previous = history[-1]
                point = (previous[0] + length * cos(angle), previous[1] + length * sin(angle))
                history.append(point)

                self._time('appened point')

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

                        # resurrect it
                        live = True

                # check for hitting detector
                if previous[0] <= self.screen <= point[0]:

                    # kill it but keep it
                    live = False
                    keep = True

                self._time('checked_status')

            # add history to self if keeping
            if keep:

                # append to self
                self.append(history)

                # print status
                if len(self) % 10 == 0:

                    # print status
                    print('{} electrons'.format(len(self), self.electrons))

                # save
                if len(self) % 10000 == 0:

                    # save
                    deposit = 'double_slit/histories_{}.json'.format(str(len(self)))
                    self._dump(self[-1000:], 'histories.json')

                    # summarize run
                    final = time()
                    percent = round(100 * self.electrons / count, 2)
                    duration = round((final - start) / 60, 2)
                    rate = round(self.electrons / duration, 0)
                    average = round(numpy.average([len(history) for history in self]), 0)
                    deviation = round(numpy.std([len(history) for history in self]), 0)

                    # print results
                    print('{} successful detections out of {} total, or {}%'.format(self.electrons, count, percent))
                    print('took {} minutes, or {} electrons / minute'.format(duration, rate))
                    print('average steps: {}, with a standard deviation of {}'.format(average, deviation))

        # summarize run
        final = time()
        percent = round(100 * self.electrons / count, 2)
        duration = round((final - start) / 60, 2)
        rate = round(self.electrons / duration, 0)
        average = round(numpy.average([len(history) for history in self]), 0)
        deviation = round(numpy.std([len(history) for history in self]), 0)

        # print results
        print('{} successful detections out of {} total, or {}%'.format(self.electrons, count, percent))
        print('took {} minutes, or {} electrons / minute'.format(duration, rate))
        print('average steps: {}, with a standard deviation of {}'.format(average, deviation))

        # save the histories
        self._dump(self, 'double_slit/histories.json')

        return None

    def populate(self):
        """Populate with saved histories rather than generating new ones.

        Arguments:
            None

        Returns:
            None

        Populates:
            self
        """

        # depopulate
        while len(self) > 1:

            # discard
            discard = self.pop()

        # open up histories
        histories = self._load('double_slit/histories.json')

        # repopulate
        [self.append(history) for history in histories]

        return None

    def spray(self, block=10000):
        """Spray electrons toward the cathode in large blocks at a time.

        Arguments:
            block=10000: block size

        Returns:
            None

        Populates:
            self
        """

        # repeat until full
        full = False
        while not full:

            # begin a block of electrons at the source
            horizontals = numpy.full(block, self.source[0]).reshape(-1, 1)
            verticals = numpy.full(block, self.source[1]).reshape(-1, 1)
            electrons = numpy.concatenate([horizontals, verticals], axis=1).reshape(block, 1, 2)

            # whittle down electrons until they hit the detector or a wall
            while len(electrons) > 0:

                # count number of surviving electrons
                survivors = len(electrons)

                # create set of random angles
                angles = random.rand(survivors) * 2 * pi

                # create set of lengths all at 1
                lengths = numpy.full(survivors, 1.0)

                # create set of random quantiles
                quantiles = random.rand(survivors)





            # set to full
            full = True

        return electrons

    def tabulate(self, precision=None, tolerance=1e-14):
        """Tabulate values of the quantile integral to be looked up during random walk generation.

        Arguments:
            precision=3: int, precision of table
            tolerance=1e-16: float, tolerance for newton-rhapson convergence

        Returns:
            None
        """

        # set precision
        precision = precision or self.precision

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

        # store lengths indexed by 10^precision * quantile
        self._dump(lengths, 'table.json')

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

    def view(self, number=1000, resolution=None):
        """View the histories.

        Arguments:
            number: int, max number of histories to plot
            resolution: int, number of histogram boxes

        Returns:
            None
        """

        # begin plot
        pyplot.clf()

        # set resolution
        resolution = resolution or self.resolution

        # plot the histograph
        self._graph()

        # set trajectory colors
        colors = ['r-', 'b-', 'c-', 'm-']
        highlights = ['w-']

        # create trajectories for the subset
        trajectories = []
        for index, history in enumerate(self[:number]):

            # use rotating colors for the majority
            color = colors[index % 4]
            if index < 1:

                # but highlights for the first few
                color = highlights[index]

            # use an exponential for quickly decaying opacity
            opacity = min([1.0, exp(-0.3 * index) + 0.01])

            # add to the trajectory
            trajectory = (history, color, opacity)
            trajectories.append(trajectory)

        # reverse the trajectories to plot highlights last
        trajectories.reverse()

        # trace the majority of trajectories
        self._trace(trajectories[:-20])

        # build the apparatus
        self._build()

        # trace the remaining trajectories
        self._trace(trajectories[-20:])

        # ignite the source
        self._ignite()

        # remove ticks
        pyplot.gca().set_xticks([])
        pyplot.gca().set_yticks([])

        # save
        pyplot.savefig('histories.png')

        return

# create instance
historian = Historian(10)
# historian.emit()
# historian.spray()
# historian.populate()
# historian.view(resolution=1000)