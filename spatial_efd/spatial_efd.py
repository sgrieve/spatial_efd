import numpy as np
import matplotlib.pyplot as plt
import shapefile as sf


def RotateContour(X, Y, rotation, centroid):
    '''
    X and Y are lists of coordinates for a closed contour taken from a shapefile
    '''

    rxs = []
    rys = []

    for nx, ny in zip(X, Y):
        rx, ry = rotatePoint(centroid, (nx, ny), rotation + 180.)
        rxs.append(rx)
        rys.append(ry)

    return rxs, rys


def NormContour(X, Y, rawCentroid):
    '''
    X and Y are the rotated contour coordinates, ready to be normed.
    '''

    # find longest axis of rotated shape
    xwidth, ywidth, xmin, ymin = getBBox(X, Y)
    if (xwidth > ywidth):
        normshape = xwidth
    elif (ywidth >= xwidth):
        normshape = ywidth

    norm_x = [(value - xmin) / normshape for value in X]
    norm_y = [(value - ymin) / normshape for value in Y]

    centroid = ((rawCentroid[0] - xmin) / normshape,
                (rawCentroid[1] - ymin) / normshape)

    return norm_x, norm_y, centroid


def ContourArea(X, Y):
    '''
    Computes the area of a closed contour
    http://paulbourke.net/geometry/polygonmesh/
    '''
    X = X + [X[0]]
    Y = Y + [Y[0]]

    Sum = 0.

    for i in xrange(len(X) - 1):
        Sum += (X[i] * Y[i + 1]) - (X[i + 1] * Y[i])

    return 0.5 * Sum


def ContourCentroid(X, Y):
    '''
    Computes the centroid of a closed contour
    http://paulbourke.net/geometry/polygonmesh/
    '''
    Area = ContourArea(X, Y)

    X = X + [X[0]]
    Y = Y + [Y[0]]

    Cx = 0.
    Cy = 0.

    for i in xrange(len(X) - 1):
        const = (X[i] * Y[i + 1]) - (X[i + 1] * Y[i])

        Cx += (X[i] + X[i + 1]) * const
        Cy += (Y[i] + Y[i + 1]) * const

    AreaFactor = (1. / (6. * Area))

    Cx *= AreaFactor
    Cy *= AreaFactor

    return (Cx, Cy)


def CalculateEFD(X, Y, harmonics=10):
    '''
    X and Y are lists of the coordinates of the normed and rotated contour.
    harmonics is the number of harmonics to be computed.

    code adapted from pyefd module.
    '''

    # pyefd takes data as y,x so build the contour like that.
    contour = np.array([(x, y) for x, y in zip(Y, X)])

    dxy = np.diff(contour, axis=0)
    dt = np.sqrt((dxy ** 2.).sum(axis=1))
    t = np.concatenate([([0., ]), np.cumsum(dt)])
    T = t[-1]

    phi = (2. * np.pi * t) / T

    coeffs = np.zeros((harmonics, 4))
    for n in xrange(1, harmonics + 1):
        const = T / (2. * n * n * np.pi * np.pi)
        phi_n = phi * n
        d_cos_phi_n = np.cos(phi_n[1:]) - np.cos(phi_n[:-1])
        d_sin_phi_n = np.sin(phi_n[1:]) - np.sin(phi_n[:-1])
        a_n = const * np.sum((dxy[:, 0] / dt) * d_cos_phi_n)
        b_n = const * np.sum((dxy[:, 0] / dt) * d_sin_phi_n)
        c_n = const * np.sum((dxy[:, 1] / dt) * d_cos_phi_n)
        d_n = const * np.sum((dxy[:, 1] / dt) * d_sin_phi_n)
        coeffs[n - 1, :] = a_n, b_n, c_n, d_n

    return coeffs


def inverse_transform(coeffs, locus=(0., 0.), n=300, harmonic=10):
    '''
    Code modified from pyefd plotting routines to generate x,y coordinates from
    forier coefficients.
    '''

    N = coeffs.shape[0]
    N_half = int(np.ceil(N / 2.))
    n_rows = 2

    t = np.linspace(0, 1., n)
    xt = np.ones((n,)) * locus[1]
    yt = np.ones((n,)) * locus[0]

    for n in xrange(harmonic):

        xt += ((coeffs[n, 0] * np.cos(2. * (n + 1.) * np.pi * t)) +
               (coeffs[n, 1] * np.sin(2. * (n + 1.) * np.pi * t)))

        yt += ((coeffs[n, 2] * np.cos(2. * (n + 1.) * np.pi * t)) +
               (coeffs[n, 3] * np.sin(2. * (n + 1.) * np.pi * t)))

        if n == harmonic - 1:
            return yt, xt


def InitPlot():
    '''
    Set up axes for a plot.
    '''
    ax = plt.gca()
    ax.axis('equal')

    return ax


def PlotEllipse(ax, x, y, color='k', width=1.):
    '''
    Plots an ellipse represented as a series of x and y coordinates on a given
    axis, ax.
    '''
    ax.plot(x, y, color, linewidth=width)


def SavePlot(ax, harmonic):
    '''
    Wrapper around the savefig method, to add a title to an axis, ax, and save
    the plot to a file.
    '''
    ax.set_title('Harmonic: {0}'.format(harmonic))
    plt.savefig('Code/test/EFD_{0}.png'.format(harmonic))
    plt.clf()


def PlotContour(ax, contour):
    '''
    Plot a contour generated from a loaded shapefile.
    '''
    ax.plot(contour[:, 1], contour[:, 0], 'b-', linewidth=2)


def AverageCoefficients(coeffList, nHarmonics):
    '''
    Average the coefficients contained in the list of coefficient arrays,
    coeffList.

    2-D particle shape averaging and comparison using Fourier descriptors:
    Powder Technology Volume 104, Issue 2, 1 September 1999, Pages 180-189
    '''
    coeffsum = np.zeros((nHarmonics, 4))
    for coeff in coeffList:
        coeffsum += coeff

    coeffsum /= float(len(coeffList))

    return coeffsum


def AverageSD(coeffList, avgcoeffs, nHarmonics):
    '''
    Use the coefficients contained in the list of coefficient arrays,
    coeffList, and the average coefficient values to compute the standard
    deviation of series of ellipses.

    2-D particle shape averaging and comparison using Fourier descriptors:
    Powder Technology Volume 104, Issue 2, 1 September 1999, Pages 180-189
    '''
    coeffsum = np.zeros((nHarmonics, 4))
    for coeff in coeffList:
        coeffsum += (coeff ** 2.)

    return (coeffsum / float(len(coeffList) - 1)) - (avgcoeffs ** 2.)


def Nyquist(contour):
    '''
    Returns the maximum number of harmonics that can be computed for a given
    contour, the nyquist freqency. Returns an integer.

    C. Costa et al. / Postharvest Biology and Technology 54 (2009) 38-47
    '''
    return contour.shape[0] / 2


def FourierPower(coeffs, contour):
    '''
    C. Costa et al. / Postharvest Biology and Technology 54 (2009) 38-47
    '''
    nyquist = Nyquist(contour)

    totalPower = 0.
    currentPower = 0.

    for n in xrange(nyquist):
            totalPower += ((coeffs[n, 0] ** 2.) + (coeffs[n, 1] ** 2.) +
                           (coeffs[n, 2] ** 2.) + (coeffs[n, 3] ** 2.)) / 2.

    for i in xrange(nyquist):
        currentPower += ((coeffs[i, 0] ** 2.) + (coeffs[i, 1] ** 2.) +
                         (coeffs[i, 2] ** 2.) + (coeffs[i, 3] ** 2.)) / 2.

        if (currentPower / totalPower) > 0.9999:
            return i + 1


def normalize_efd(coeffs, size_invariant=True):
    '''
    Taken and modified to run in python 2.7 from pyefd.
    '''
    # Make the coefficients have a zero phase shift from
    # the first major axis. Theta_1 is that shift angle.
    theta_1 = (0.5 * np.arctan2(2. * ((coeffs[0, 0] * coeffs[0, 1]) +
               (coeffs[0, 2] * coeffs[0, 3])),
              ((coeffs[0, 0] ** 2.) -
               (coeffs[0, 1] ** 2.) +
               (coeffs[0, 2] ** 2.) -
               (coeffs[0, 3] ** 2.))))

    # Rotate all coefficients by theta_1.
    for n in xrange(1, coeffs.shape[0] + 1):
        coeffs[n - 1, :] = np.dot(np.array([[coeffs[n - 1, 0],
                                  coeffs[n - 1, 1]], [coeffs[n - 1, 2],
                                                      coeffs[n - 1, 3]]]),
                                  np.array([[np.cos(n * theta_1),
                                           -np.sin(n * theta_1)],
                                           [np.sin(n * theta_1),
                                           np.cos(n * theta_1)]])).flatten()

    # Make the coefficients rotation invariant by rotating so that
    # the semi-major axis is parallel to the x-axis.
    psi_1 = np.arctan2(coeffs[0, 2], coeffs[0, 0])
    psi_r = np.array([[np.cos(psi_1), np.sin(psi_1)],
                     [-np.sin(psi_1), np.cos(psi_1)]])

    # Rotate all coefficients by -psi_1.
    for n in xrange(1, coeffs.shape[0] + 1):
        coeffs[n - 1, :] = psi_r.dot(np.array([[coeffs[n - 1, 0],
                                     coeffs[n - 1, 1]], [coeffs[n - 1, 2],
                                                         coeffs[n - 1,
                                                                3]]])).flatten()

    if size_invariant:
        # Obtain size-invariance by normalizing.
        coeffs /= np.abs(coeffs[0, 0])

    return coeffs


def calculate_dc_coefficients(contour):
    '''
    Contour is a list of Y, X coordinate tuples.
    code adapted from pyefd module.
    '''
    dxy = np.diff(contour, axis=0)
    dt = np.sqrt((dxy ** 2.).sum(axis=1))
    t = np.concatenate([([0., ]), np.cumsum(dt)])
    T = t[-1]

    diff = np.diff(t ** 2.)
    xi = np.cumsum(dxy[:, 0]) - (dxy[:, 0] / dt) * t[1:]
    A0 = (1. / T) * np.sum(((dxy[:, 0] / (2. * dt)) * diff) + xi * dt)
    delta = np.cumsum(dxy[:, 1]) - (dxy[:, 1] / dt) * t[1:]
    C0 = (1. / T) * np.sum(((dxy[:, 1] / (2. * dt)) * diff) + delta * dt)

    # A0 and CO relate to the first point of the contour array as origin.
    # Adding those values to the coefficients to make them relate to true origin
    return (contour[0, 1] + C0, contour[0, 0] + A0)


def LoadGeometries(filename):
    '''
    Takes a filename and uses shapefile to load it, returning a list of
    shapefile.ShapeRecord instances.

    this list can be iterated over, passing the shape individual shape instances
    to ProcessGeometry() one by one.
    '''

    return sf.shapeRecords()


def ProcessGeometry(shape):
    '''
    Method which takes a single shape instance from a shapefile
    eg shp.Reader('shapefile.shp').shapeRecords()[n]
    where n is the index of the shape within a multipart geometry, and returns
    a list of x coordinates, a list of y coordinates, contour (a list of
    [x,y] coordinate pairs, normalized about the shape's centroid) and
    the normalized coordinate centroid.

    This handles all the geometry conversion that may be needed by the rest of
    the EFD code.
    '''
    x = []
    y = []

    for point in shape.shape.points:
        x.append(point[0])
        y.append(point[1])

    centroid = ContourCentroid(x, y)
    X, Y, NormCentroid = NormContour(x, y, centroid)

    contour = np.array([(a, b) for a, b in zip(Y, X)])

    return x, y, contour, NormCentroid


def rotatePoint(centerPoint, point, angle):
    """
    Rotates a point around another centerPoint. Angle is in degrees.
    Rotation is counter-clockwise
    https://gist.github.com/somada141/d81a05f172bb2df26a2c
    """
    angle = math.radians(angle)
    temp_point = point[0] - centerPoint[0], point[1] - centerPoint[1]
    temp_point = (temp_point[0] * math.cos(angle) - temp_point[1] *
                  math.sin(angle), temp_point[0] * math.sin(angle) +
                  temp_point[1] * math.cos(angle))

    temp_point = temp_point[0] + centerPoint[0], temp_point[1] + centerPoint[1]
    return temp_point[0], temp_point[1]


def getBBox(x, y):
    xmin = min(x)
    ymin = min(y)

    xmax = max(x)
    ymax = max(y)

    return xmax - xmin, ymax - ymin, xmin, ymin
