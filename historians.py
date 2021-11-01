# historians.py for the historian class to keep track of particle histories

# import system tools
import os, sys
import json

# import time
from time import time, sleep

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

    def __init__(self, directory='testo', electrons=10, wave=1000, mode='sin2', status=True, statusii=True):
        """Initialize a historian instance.

        Arguments:
            directory: str, name of directory off root to place files
            electrons: int, number of successful electrons
            wave: int, number of electrons per file
            mode: str: distribution mode
            status: boolean, first slit open?
            statusii: boolean, second slit open?

        Returns:
            None
        """

        # set current directory
        self.directory = directory
        self._make(directory)

        # current time
        self.now = time()

        # set number of total electrons to capture
        self.electrons = int(electrons)

        # set number of electrons per wave
        self.wave = int(wave)

        # set up apparatus
        self.top = 40
        self.bottom = -40
        self.back = -20
        self.source = (0, 0)
        self.divider = 40
        self.screen = 100

        # set clarity
        self.clarity = 0.2

        # configure slits
        self.gap = 0.5
        self.space = 25
        self.statuses = [bool(status), bool(statusii)]
        self.slits = []
        self._configure()

        # use weighting at detector screen?
        self.weightings = True

        # set histogram resolution
        self.resolution = 10000

        # set tabulation precision
        self.precision = 3

        # define probability distribution
        self.mode = mode
        self.distribution = None
        self.quantile = None
        self.zero = None
        self.slope = None
        self.guess = None
        self.coda = None
        self.inverse = None
        self.zeroii = None
        self.slopeii = None
        self._formulate()

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

    def _crank(self, quantile, zero, slope, guess=1.0, trials=5, tolerance=1e-14):
        """Crank through Newton Rhapson to solve the quantile equation for a segment length.

        Arguments:
            quantile: float, the quantile to pinpoint
            zero: function object
            slope: function object
            guess=1.0: float, initial guess
            trials: int, number guesses
            tolerance=1e-12: tolerance for closesness to zero

        Returns:
            float, better guess for segment length
        """

        # add random number to guess for several starting guesses
        guesses = [guess + (random() - 0.5) * 0.01 for _ in range(trials)]

        # get roots
        roots = []

        # for each guess
        for guess in guesses:

            # evaluate the function
            evaluation = zero(guess, quantile)
            derivative = slope(guess, quantile)

            # while the function evaluates to outside the tolerance
            while abs(evaluation) > tolerance:

                # adjust using newton's formula
                guess = guess - evaluation / derivative

                # evaluate the function
                evaluation = zero(guess, quantile)
                derivative = slope(guess, quantile)

            # add to roots
            roots.append(guess)

        # whittle down to real numbers += 10 (should cover to 3pi)
        roots = [root for root in roots if numpy.isfinite(root)]
        roots = [root for root in roots if abs(root) < 10]

        # average remaining for final value
        root = numpy.average(roots)

        return root

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

    def _formulate(self):
        """Define distribution functions.

        Arguments:
            mode: str,  particular mode

        Returns:
            None

        Populates:
            self.distribution
            self.quantile
            self.zero
            self.slope
        """

        # define sin^2 distribution functions
        if self.mode == 'sin2':

            # define probability distribution functions
            self.distribution = lambda t:  2 * cos(pi * t) ** 2
            self.quantile = lambda t: t - 0.5 + sin(2 * pi * t) / (2 * pi)

            # prepare functions for newton-rhapson
            self.zero = lambda t, q: t - q - 0.5 + sin(2 * pi * t) / (2 * pi)
            self.slope = lambda t, q: 1 + cos(2 * pi * t)

            # add coda
            self.coda = lambda t: t

            # add newtton rhpapson parameters
            self.zeroii = lambda t, x: t - x
            self.slopeii = lambda t, x: 1

            # set guess
            self.guess = 1.0

        # define cycloid distribution functions
        if self.mode == 'cycloid':

            # define probability distribution functions
            self.distribution = lambda t:  (25 / (4 * pi)) * (1 + cos(t))
            self.quantile = lambda t: (1 / pi) * ((t - pi) / 2 + (-sin(2 * t) / 4))

            # prepare functions for newton-rhapson
            self.zero = lambda t, q: (1 / pi) * ((t - pi) / 2 + (-sin(2 * t) / 4)) - q
            self.slope = lambda t, q: (1 / pi) * ((1 / 2) + (-cos(2 * t) / 2))

            # add coda
            self.coda = lambda t: (4 / 25) * (t - sin(t))

            # add coda newton parameterrs
            self.zeroii = lambda t, x: (4 / 25) * (t - sin(t)) - x
            self.slopeii = lambda t, x: (4 / 25) * (1 - cos(t))

            # set guess
            self.guess = 1.0

        return None

    def _graph(self, resolution):
        """Plot the histograph.

        Arguments:
            resolution: int, number of bins

        Returns:
            None
        """

        # make bins
        self._stamp('making bins...')
        chunk = (self.top - self.bottom) / resolution
        bins = [(self.bottom + index * chunk, self.bottom + (index + 1) * chunk) for index in range(resolution)]

        # populate bins according to weighting scheme (-2 before (?))
        self._stamp('populating bins...')
        population = [sum([self._weigh(member) for member in self if bin[0] <= member[-1][1] < bin[1]]) for bin in bins]

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

    def _make(self, folder):
        """Make a folder in the directory if not yet made.

        Arguments:
            folder: str, directory path

        Returns:
            None
        """

        # try to
        try:

            # create the directory
            os.mkdir(folder)

        # unless directory already exists, or none was given
        except (FileExistsError, FileNotFoundError, PermissionError):

            # in which case, nevermind
            pass

        return None

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

    def _weigh(self, member):

        # set default weight to the destribution at length 1
        weight = self.distribution(1.0)

        # if the leg is non classical (not strictly horizontal) and weightings is set
        if self.weightings and member[-1][1] != member[-2][1]:

            # measure the length of the leg and apply distribution
            length = self._measure(member[-1], member[-2])
            coda = self._crank(length, self.zeroii, self.slopeii)
            weight = self.distribution(coda)

        return weight

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

        # convert to parameterization
        parameterization = [self._crank(horizontal, self.zeroii, self.slopeii) for horizontal in horizontals]

        # calculate functions
        distributions = numpy.array([self.distribution(horizontal) for horizontal in parameterization])
        quantiles = numpy.array([self.quantile(horizontal) for horizontal in parameterization])

        # set up plot
        pyplot.clf()

        # plot functions
        pyplot.plot(horizontals, distributions, 'b--')
        pyplot.plot(horizontals, quantiles, 'g--')

        # save figure
        pyplot.savefig('verifications/distribution_{}.png'.format(self.mode))

        return None

    def populate(self, number=1000, entirety=5):
        """Populate with saved histories rather than generating new ones.

        Arguments:
            number: int, number of files
            entirety: int, number of files to grab in total, otherwise just last steps

        Returns:
            None

        Populates:
            self
        """

        # depopulate
        self._stamp('populating...')
        while len(self) > 1:

            # discard
            self.pop()

        # get all histories
        histories = []
        waves = os.listdir(self.directory)
        waves.sort(reverse=True)
        waves = [wave for wave in waves if '.png' not in wave][:number]
        for index, wave in enumerate(waves):

            # open up histories
            path = '{}/{}'.format(self.directory, wave)
            quiver = self._load(path)

            # if beyond entirety:
            if index > entirety:

                # subset to last few paths
                quiver = [arrow[-10:] for arrow in quiver]

            # add to histories
            histories += quiver

        # repopulate
        [self.append(history) for history in histories]

        # print length
        print('{} trajectories.'.format(len(histories)))

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

        # open up tabulated values and convert to dict
        table = self._load('tables/table_{}.json'.format(self.mode))
        table = {index: entry for index, entry in enumerate(table)}

        # begin the clock and count
        start = time()
        count = 0
        successes = 0
        accumulation = 0

        # scan all files in directory
        waves = os.listdir(self.directory)
        total = sum([int(wave.split('_')[2].split('.')[0]) for wave in waves if 'png' not in wave])
        wave = len(waves)

        # repeat until full
        while total < self.electrons:

            # increase count
            count += 1

            # begin a block of electrons at the source
            horizontals = numpy.full(block, self.source[0]).reshape(-1, 1)
            verticals = numpy.full(block, self.source[1]).reshape(-1, 1)
            electrons = numpy.concatenate([horizontals, verticals], axis=1).reshape(-1, 1, 2)

            # propagate electrons until they hit the detector or a wall
            while len(electrons) > 0:

                # count number of surviving electrons
                survivors = len(electrons)

                # create set of random lengths, with maximum set by clarity parameter
                lengths = rand(survivors) * self.clarity

                # create classical steps toward screen
                horizontals = numpy.add(electrons[:, -1, 0], lengths).reshape(-1, 1)
                verticals = electrons[:, -1, 1].reshape(-1, 1)
                classicals = numpy.concatenate([horizontals, verticals], axis=1).reshape(-1, 1, 2)

                # add to electron trajectories
                electrons = numpy.concatenate([electrons, classicals], axis=1)

                # keep all the electrons that hit the detector but remove from live electrons
                detections = electrons[electrons[:, -1, 0] > self.screen]
                electrons = electrons[electrons[:, -1, 0] < self.screen]

                # add to successes
                successes += len(detections)
                accumulation += len(detections)

                # cut detections off at the screen
                if len(detections) > 0:

                    # determine cutoff point and append
                    detections[:, -1, 0] = self.screen
                    [self.append(history) for history in detections]

                # find those that span the divider and remove them
                spanning = (electrons[:, -2, 0] < self.divider) != (electrons[:, -1, 0] < self.divider)
                electrons = electrons[~spanning]

                # recount number of surviving electrons
                survivors = len(electrons)

                # create set of random angles for qantum paths
                angles = rand(survivors) * 2 * pi

                # create set of random ints
                indices = randint(1001, size=(survivors,))

                # get lengths from table
                lengths = numpy.array([table[index] for index in indices])

                # make new random walk
                horizontals = numpy.add(electrons[:, -1, 0], numpy.multiply(lengths, cos(angles))).reshape(-1, 1)
                verticals = numpy.add(electrons[:, -1, 1], numpy.multiply(lengths, sin(angles))).reshape(-1, 1)
                walk = numpy.concatenate([horizontals, verticals], axis=1).reshape(-1, 1, 2)

                # add to electron trajectories
                electrons = numpy.concatenate([electrons, walk], axis=1)

                # prune off those that hit the top and bottom
                electrons = electrons[electrons[:, -1, 1] < self.top]
                electrons = electrons[electrons[:, -1, 1] > self.bottom]

                # prune off those that hit the back
                electrons = electrons[electrons[:, -1, 0] > self.back]

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

                # find those that span the divider and remove them
                spanning = (electrons[:, -2, 0] < self.divider) != (electrons[:, -1, 0] < self.divider)
                spanners = electrons[spanning]
                electrons = electrons[~spanning]

                # find the slopes from the approaching points to all slit edges
                approaches = []
                edges = [entry for slit in self.slits for entry in slit]
                for edge in edges:

                    # create the array
                    approach = (edge - spanners[:, -2, 1]) / (self.divider - spanners[:, -2, 0])
                    approaches.append(approach)

                # find the slopes from the departing points to all slit edges
                departures = []
                for edge in edges:

                    # create the array
                    departure = (spanners[:, -1, 1] - edge) / (spanners[:, -1, 0] - self.divider)
                    departures.append(departure)

                # calculate all differences between departing and approaching slopes (curvatures)
                curvatures = [numpy.subtract(approach, departure) for approach, departure in zip(approaches, departures)]

                # if the product of all curvatures is negative, the electron went through a slit
                product = numpy.prod(numpy.vstack(curvatures), axis=0)
                passers = spanners[product < 0]

                # add passers back to electrons
                electrons = numpy.concatenate([electrons, passers], axis=0)

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
                deposit = '{}/histories_{}_{}.json'.format(self.directory, str(wave).zfill(4), successes)
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
                print('calculating quantile {}, {} of {}'.format(quantile, index, len(quantiles)))

            # calculate length and make a string
            length = self._crank(quantile, self.zero, self.slope, guess=self.guess, tolerance=tolerance)
            length = self.coda(length)
            lengths.append(length)

        # store lengths indexed by 10^precision * quantile
        self._dump(lengths, 'tables/table_{}.json'.format(self.mode))

        # print table
        pyplot.clf()

        # plot table
        indices = [index for index, _ in enumerate(lengths)]
        pyplot.plot(indices, lengths, 'g--')
        pyplot.savefig('verifications/table_{}'.format(self.mode))
        pyplot.clf()

        return None

    def verify(self, trials=10000):
        """Test the distribution resulting from the inverse quantile function and lookup table.

        Arguments:
            None

        Returns:
            None
        """

        # load in table
        table = self._load('tables/table_{}.json'.format(self.mode))

        # find all lengths based on random number selection
        lengths = [table[int(random() * 1000)] for _ in range(trials)]

        # bin them to draw a histogram
        middles, heights = self._bin(lengths, resolution=100, start=0.5, finish=1.5)

        # calculate distribution function
        chunk = 0.01
        horizontals = [number * chunk + 0.5 for number in range(101)]

        # convert to parameterization
        parameterization = [self._crank(horizontal, self.zeroii, self.slopeii) for horizontal in horizontals]

        # calculate distribution
        distributions = numpy.array([self.distribution(horizontal) for horizontal in parameterization])

        # start plot
        pyplot.clf()

        # plot bars
        for middle, height in zip(middles, heights):

            # plot a bar
            pyplot.plot([middle, middle], [0, height], 'g-', linewidth=1)

        # plot distribution
        pyplot.plot(horizontals, distributions, 'b--')

        # save fig
        pyplot.savefig('verifications/verification_{}.png'.format(self.mode))

        return None

    def view(self, number=5000, resolution=128):
        """View the histories.

        Arguments:
            number: int, max number of histories to plot
            resolution: int, number of histogram boxes

        Returns:
            None
        """

        # begin plot
        self._stamp('creating plot...')
        pyplot.clf()

        # set resolution
        resolution = resolution or self.resolution

        # plot the histograph
        self._stamp('generating histogram...')
        self._graph(resolution)

        # set trajectory colors
        colors = ['r-', 'b-', 'c-', 'm-']
        highlights = ['w-']

        # create trajectories for the subset
        trajectories = []
        self._stamp('gathering trajectories...')
        for index, history in enumerate(self[:number]):

            # use rotating colors for the majority
            color = colors[index % 4]
            if index < 1:

                # but highlights for the first few
                color = highlights[index]

            # use an exponential for quickly decaying opacity
            opacity = min([1.0, exp(-0.1 * index) + 0.01])

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

        print(dir(pyplot.gca().yaxis))
        print(pyplot.gca().yaxis.get_scale())

        # save
        deposit = 'experiments/{}_{}_{}.png'.format(self.directory.split('/')[-1], len(self), resolution)
        pyplot.savefig(deposit)

        self._stamp('{} saved.'.format(deposit))

        return


# grab potential arguments
arguments = [argument for argument in sys.argv[1:]]

# if there ae arguments
if len(arguments) > 0:

    # create instance
    historian = Historian(*arguments)

    # spray until enough electrons hit detector
    historian.spray()

    # gather all histories
    historian.populate()

    # tally histogram and plot
    historian.view()