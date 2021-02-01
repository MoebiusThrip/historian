# historians.py for the historian class to keep track of particle histories

# import reload
from importlib import reload

# import system tools
import os, sys
import json

# import time
from time import time, sleep

# import random

# import trig functions
from math import sin, cos, pi, exp, log

# import numpy
import numpy
from numpy.random import random, rand, randint
from numpy import sin, cos, pi, exp, log, sqrt

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

    def __init__(self, directory, electrons=100, wave=10000, status=True, statusii=True):
        """Initialize a historian instance.

        Arguments:
            electrons: int, number of successful electrons
            status: boolean, first slit open?
            statusii: boolean, second slit open?

        Returns:
            None
        """

        # set current directory
        self.directory = directory

        # current time
        self.now = time()

        # set number of total electrons to capture
        self.electrons = electrons

        # set number of electrons per wave
        self.wave = wave

        # set up apparatus
        self.top = 20
        self.bottom = -20
        self.back = -20
        self.divider = 0
        self.screen = 20
        self.source = (-10, 0)

        # configure slits
        self.gap = 0.5
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
        height = point[0][1] * (pointii[0][0] - horizontal) + pointii[0][1] * (horizontal - point[0][0]) / (pointii[0][0] - point[0][0])

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

        # make bins
        self._stamp('making bins...')
        chunk = (self.top - self.bottom) / self.resolution
        bins = [(self.bottom + index * chunk, self.bottom + (index + 1) * chunk) for index in range(self.resolution)]

        # populate bins
        self._stamp('populating bins...')
        population = [sum([self.distribution(self._measure(member[-1], member[-2])) for member in self if bin[0] <= member[-2][1] < bin[1]]) for bin in bins]
        #population = [len([member for member in self if bin[0] < self._cross(self.screen, member[0], member[1]) < bin[1]]) for bin in bins]

        # normalize population
        self._stamp('normalizing...')
        maximum = max(population)
        population = [entry * 10 / maximum for entry in population]

        # plot it
        self._stamp('plotting...')
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

    def _measure(self, point, pointii):
        """Measure the distance between two points.

        Arguments:
            point: tuple of floats
            pointii: tuple of floats

        Returns:
            float
        """

        # calculate length
        length = sqrt((pointii[0] - point[0]) ** 2 + (pointii[1] - point[1]) ** 2)

        return length

    def _stamp(self, message):
        """Start timing a block of code, and print results with a message.

        Arguments:
            message: str

        Returns:
            None
        """

        # get final time
        final = time()

        # calculate duration and reset time
        duration = final - self.now
        self.now = final

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

            self._stamp('starting...')

            # random walk until it hits something
            live = True
            keep = False
            while live:

                self._stamp('starting electron...')

                # generate more random numbers
                if len(randoms) < 10:

                    # generate more random numbers
                    randoms = random(size=1000000).tolist() + randoms

                # pick angle at random
                angle = randoms.pop() * 2 * pi

                # look up length from the table using a randomly generated quantile
                index = int(randoms.pop() * 10 ** self.precision)
                length = table[index]

                self._stamp('made randoms')

                # create new point
                previous = history[-1]
                point = (previous[0] + length * cos(angle), previous[1] + length * sin(angle))
                history.append(point)

                self._stamp('appened point')

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

                self._stamp('checked_status')

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

        self._stamp('populating...')

        # depopulate
        while len(self) > 1:

            # discard
            discard = self.pop()

        # get all histories
        histories = []
        waves = os.listdir(self.directory)
        waves = [wave for wave in waves if '.png' not in wave]
        for wave in waves:

            # open up histories
            path = '{}/{}'.format(self.directory, wave)
            histories += self._load(path)

        # repopulate
        [self.append(history) for history in histories]

        # print length
        print('{} trajectories.'.format(len(histories)))

        return None

    def spatter(self, block=1000, steps=1000, tolerance=1e-14):
        """Spatter electrons toward the cathode in large blocks at a time.

        Arguments:
            block=1000: block size
            steps=1000: number of steps in path

        Returns:
            None

        Populates:
            self
        """

        # open up tabulated quartile values
        table = self._load('table.json')

        # begin the clock and count
        start = time()
        count = 0
        successes = 0
        accumulation = 0

        # scan all files in directory
        waves = os.listdir(self.directory)
        total = sum([int(wave.split('_')[1]) for wave in waves if 'png' not in wave])
        wave = len(waves)

        # self._stamp('begin collection...')

        # repeat until full
        while total < self.electrons:

            # print
            # self._stamp('spraying {} electrons...'.format(block))

            # increase count
            count += 1

            # # begin a block of electrons at the source
            # horizontals = numpy.full(block, self.source[0]).reshape(-1, 1)
            # verticals = numpy.full(block, self.source[1]).reshape(-1, 1)
            # electrons = numpy.concatenate([horizontals, verticals], axis=1).reshape(-1, 1, 2)

            # begin a block of electrons
            electrons = numpy.zeros((block, steps, 2))

            # print(electrons.shape)

            # begin each electron at the source
            step = 0
            electrons[:, 0] = self.source

            # self._stamp('began block')

            # propagate electrons until they hit the detector or a wall
            while len(electrons) > 0:

                # print(step, len(electrons))

                # increase step
                step += 1

                # check for last step
                if step == len(electrons[0]) - 1:

                    # add new block
                    news = numpy.zeros((len(electrons), steps, 2))
                    electrons = numpy.concatenate([electrons, news], axis=1)

                # # increase stone
                # stones += 1
                # if stones % 400 == 0:
                #
                #     # print status
                #     print('{} steps so far, {} electrons left...'.format(stones, len(electrons)))

                # count number of surviving electrons
                survivors = len(electrons)

                # create set of random angles
                angles = rand(survivors) * 2 * pi
                #
                # # create set of lengths all at 1
                # lengths = numpy.full(survivors, 1.0)

                # self._stamp('made angles')

                # # create set of random quantiles
                # quantiles = rand(survivors)
                #
                # # evaluate zeros for all quantiles
                # zeros = self.zero(lengths, quantiles)
                #
                # # check for values above tolerance
                # while numpy.any(zeros > tolerance):
                #
                #     # use newton rhapson to get closer
                #     slopes = self.slope(lengths, quantiles)
                #
                #     # make new lengths
                #     lengths = lengths - zeros / slopes
                #
                #     # calculate new zeros
                #     zeros = self.zero(lengths, quantiles)

                # create set of random integers to use as indices
                indices = randint(1001, size=(survivors,))

                # self._stamp('made random ints')

                # get lengths from table
                lengths = numpy.array([table[index] for index in indices])

                # self._stamp('checked lengths')

                # make new random walk
                horizontals = numpy.add(electrons[:, step - 1, 0], numpy.multiply(lengths, cos(angles)))
                verticals = numpy.add(electrons[:, step - 1, 1], numpy.multiply(lengths, sin(angles)))
                #walk = numpy.concatenate([horizontals, verticals], axis=1)
                    #.reshape(-1, 1, 2)
                #
                # self._stamp('made new steps')

                # # add to electron trajectories
                # cube = numpy.zeros((electrons.shape[0], electrons.shape[1] + 1, 2))
                # cube[:, :len(electrons[0]), :2] = electrons
                # cube[:, len(electrons[0]):, :2] = walk
                # electrons = cube

                # add new step
                electrons[:, step, 0] = horizontals.reshape(1, -1)
                electrons[:, step, 1] = verticals.reshape(1, -1)
                # self._stamp('added to current stash')

                # prune off those that hit the top and bottom
                electrons = electrons[electrons[:, step, 1] < self.top]
                electrons = electrons[electrons[:, step, 1] > self.bottom]

                # self._stamp('pruned off top and bottom')

                # prune off those that hit the back
                electrons = electrons[electrons[:, step, 0] > self.back]

                # self._stamp('pruned off back')

                # keep all the electrons that hit the detector but remove from live electrons
                detections = electrons[electrons[:, step, 0] > self.screen]
                electrons = electrons[electrons[:, step, 0] < self.screen]

                # add to successes
                successes += len(detections)
                accumulation += len(detections)

                # cut detections off at the screen
                if len(detections) > 0:

                    # determine cutoff point and append
                    detections[:, step] = numpy.array([self.screen, self._cross(self.screen, detections[:, step - 1], detections[:, step])])
                    [self.append(history[:step + 1]) for history in detections]

                # self._stamp('getting cutoff points')

                # find those that span the divider and remove them
                spanning = (electrons[:, step - 1, 0] < self.divider) != (electrons[:, step, 0] < self.divider)
                spanners = electrons[spanning]
                electrons = electrons[~spanning]

                # self._stamp('finding spanners')

                # find the slopes from the approaching points to all slit edges
                approaches = []
                edges = [entry for slit in self.slits for entry in slit]
                for edge in edges:

                    # create the array
                    approach = (edge - spanners[:, step - 1, 1]) / (self.divider - spanners[:, step - 1, 0])
                    approaches.append(approach)

                # self._stamp('calculating approaches')

                # find the slopes from the departing points to all slit edges
                departures = []
                for edge in edges:

                    # create the array
                    departure = (spanners[:, step, 1] - edge) / (spanners[:, step, 0] - self.divider)
                    departures.append(departure)

                # self._stamp('calculating departures')

                # calculate all differences between departing and approaching slopes (curvatures)
                curvatures = [numpy.subtract(approach, departure) for approach, departure in zip(approaches, departures)]

                # self._stamp('calculating curvatures')

                # if the product of all curvatures is negative, the electron went through a slit
                product = numpy.prod(numpy.vstack(curvatures), axis=0)
                passers = spanners[product < 0]

                # self._stamp('found spanners')

                # add passers back to electrons
                electrons = numpy.concatenate([electrons, passers], axis=0)

                # self._stamp('concatenating total')

            # summarize block
            final = time()
            emissions = count * block
            percent = round(100 * accumulation / emissions, 2)
            duration = round((final - start) / 60, 2)
            rate = round(accumulation / duration, 0)
            average = round(numpy.average([len(history) for history in self]), 0)
            deviation = round(numpy.std([len(history) for history in self]), 0)

            # print results
            if count % 2 == 0:

                # pirnt
                print('{} successful detections out of {} total, or {}%'.format(accumulation, emissions, percent))
                print('took {} minutes, or {} electrons / minute'.format(duration, rate))
                print('average steps: {}, with a standard deviation of {}'.format(average, deviation))
                sleep(2)

            # save file
            if successes > self.wave:

                # get histories
                histories = []
                while len(self) > 0:

                    # depopulate
                    histories.append(self.pop())

                # save file
                deposit = '{}/histories_{}_{}.json'.format(self.directory, successes, wave)
                histories = [history.tolist() for history in histories]
                self._dump(histories, deposit)

                # increase wave
                wave += 1

                # add to total and reset count
                total += successes
                successes = 0

        return None

    def spray(self, block=1000, tolerance=1e-14):
        """Spray electrons toward the cathode in large blocks at a time.

        Arguments:
            block=10000: block size

        Returns:
            None

        Populates:
            self
        """

        # open up tabulated values
        table = self._load('table.json')

        # begin the clock and count
        start = time()
        count = 0
        successes = 0
        accumulation = 0

        # scan all files in directory
        waves = os.listdir(self.directory)
        total = sum([int(wave.split('_')[1]) for wave in waves if 'png' not in wave])
        wave = len(waves)

        # self._stamp('begin collection...')

        # repeat until full
        while total < self.electrons:

            # print
            # self._stamp('spraying {} electrons...'.format(block))

            # increase count
            count += 1

            # begin a block of electrons at the source
            horizontals = numpy.full(block, self.source[0]).reshape(-1, 1)
            verticals = numpy.full(block, self.source[1]).reshape(-1, 1)
            electrons = numpy.concatenate([horizontals, verticals], axis=1).reshape(-1, 1, 2)

            # self._stamp('began block')

            # propagate electrons until they hit the detector or a wall
            stones = 0
            while len(electrons) > 0:

                # # increase stone
                # stones += 1
                # if stones % 400 == 0:
                #
                #     # print status
                #     print('{} steps so far, {} electrons left...'.format(stones, len(electrons)))

                # count number of surviving electrons
                survivors = len(electrons)

                # create set of random angles
                angles = rand(survivors) * 2 * pi
                #
                # # create set of lengths all at 1
                # lengths = numpy.full(survivors, 1.0)

                # self._stamp('made angles')

                # # create set of random quantiles
                # quantiles = rand(survivors)
                #
                # # evaluate zeros for all quantiles
                # zeros = self.zero(lengths, quantiles)
                #
                # # check for values above tolerance
                # while numpy.any(zeros > tolerance):
                #
                #     # use newton rhapson to get closer
                #     slopes = self.slope(lengths, quantiles)
                #
                #     # make new lengths
                #     lengths = lengths - zeros / slopes
                #
                #     # calculate new zeros
                #     zeros = self.zero(lengths, quantiles)

                # create set of random ints
                indices = randint(1001, size=(survivors,))

                # self._stamp('made random ints')

                # get lengths from table
                lengths = numpy.array([table[index] for index in indices])

                # self._stamp('checked lengths')

                # make new random walk
                horizontals = numpy.add(electrons[:, -1, 0], numpy.multiply(lengths, cos(angles))).reshape(-1, 1)
                verticals = numpy.add(electrons[:, -1, 1], numpy.multiply(lengths, sin(angles))).reshape(-1, 1)
                walk = numpy.concatenate([horizontals, verticals], axis=1).reshape(-1, 1, 2)
                #
                # self._stamp('made new steps')

                # # add to electron trajectories
                # cube = numpy.zeros((electrons.shape[0], electrons.shape[1] + 1, 2))
                # cube[:, :len(electrons[0]), :2] = electrons
                # cube[:, len(electrons[0]):, :2] = walk
                # electrons = cube

                electrons = numpy.concatenate([electrons, walk], axis=1)
                # self._stamp('added to current stash')

                # prune off those that hit the top and bottom
                electrons = electrons[electrons[:, -1, 1] < self.top]
                electrons = electrons[electrons[:, -1, 1] > self.bottom]

                # self._stamp('pruned off top and bottom')

                # prune off those that hit the back
                electrons = electrons[electrons[:, -1, 0] > self.back]

                # self._stamp('pruned off back')

                # keep all the electrons that hit the detector but remove from live electrons
                detections = electrons[electrons[:, -1, 0] > self.screen]
                electrons = electrons[electrons[:, -1, 0] < self.screen]

                # add to successes
                successes += len(detections)
                accumulation += len(detections)

                # cut detections off at the screen
                if len(detections) > 0:

                    # determine cutoff point and append
                    detections[:, -1] = numpy.array([self.screen, self._cross(self.screen, detections[:, -2], detections[:, -1])])
                    [self.append(history) for history in detections]

                # self._stamp('getting cutoff points')

                # find those that span the divider and remove them
                spanning = (electrons[:, -2, 0] < self.divider) != (electrons[:, -1, 0] < self.divider)
                spanners = electrons[spanning]
                electrons = electrons[~spanning]

                # self._stamp('finding spanners')

                # find the slopes from the approaching points to all slit edges
                approaches = []
                edges = [entry for slit in self.slits for entry in slit]
                for edge in edges:

                    # create the array
                    approach = (edge - spanners[:, -2, 1]) / (self.divider - spanners[:, -2, 0])
                    approaches.append(approach)

                # self._stamp('calculating approaches')

                # find the slopes from the departing points to all slit edges
                departures = []
                for edge in edges:

                    # create the array
                    departure = (spanners[:, -1, 1] - edge) / (spanners[:, -1, 0] - self.divider)
                    departures.append(departure)

                # self._stamp('calculating departures')

                # calculate all differences between departing and approaching slopes (curvatures)
                curvatures = [numpy.subtract(approach, departure) for approach, departure in zip(approaches, departures)]

                # self._stamp('calculating curvatures')

                # if the product of all curvatures is negative, the electron went through a slit
                product = numpy.prod(numpy.vstack(curvatures), axis=0)
                passers = spanners[product < 0]

                # self._stamp('found spanners')

                # add passers back to electrons
                electrons = numpy.concatenate([electrons, passers], axis=0)

                # self._stamp('concatenating total')

            # summarize block
            final = time()
            emissions = count * block
            percent = round(100 * accumulation / emissions, 2)
            duration = round((final - start) / 60, 2)
            rate = round(accumulation / duration, 0)
            average = round(numpy.average([len(history) for history in self]), 0)
            deviation = round(numpy.std([len(history) for history in self]), 0)

            # print results
            if count % 100 == 0:

                # prnt
                print('{} successful detections out of {} total, or {}%'.format(accumulation, emissions, percent))
                print('took {} minutes, or {} electrons / minute'.format(duration, rate))
                print('average steps: {}, with a standard deviation of {}'.format(average, deviation))

            # save file
            if successes > self.wave:

                # get histories
                histories = []
                while len(self) > 0:

                    # depopulate
                    histories.append(self.pop())

                # save file
                deposit = '{}/histories_{}_{}.json'.format(self.directory, successes, wave)
                histories = [history.tolist() for history in histories]
                self._dump(histories, deposit)

                # increase wave
                wave += 1

                # add to total and reset count
                total += successes
                successes = 0

        return None

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

    def view(self, number=1000, resolution=100):
        """View the histories.

        Arguments:
            number: int, max number of histories to plot
            resolution: int, number of histogram boxes

        Returns:
            None
        """

        self._stamp('creating plot...')

        # begin plot
        pyplot.clf()

        # set resolution
        resolution = resolution or self.resolution

        self._stamp('generating histogram...')

        # plot the histograph
        self._graph()

        # set trajectory colors
        colors = ['r-', 'b-', 'c-', 'm-']
        highlights = ['w-']

        self._stamp('gathering trajectories...')

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

        self._stamp('tracing trajectories....')

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
        deposit = '{}/histories.png'.format(self.directory)
        pyplot.savefig(deposit)

        return

# create instance
historian = Historian('double', 1000000, 1000)
historian.spray()
# historian.emit()
# historian.spray()
# historian.spatter()
historian.populate()
historian.view()