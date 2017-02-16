#!/usr/bin/env python
import numpy as np
import shapefile as sf
import matplotlib.pyplot as plt


def RotateContour(X, Y, rotation, centroid):
    '''
    Rotates a contour about a point by a given amount expressed in degrees.

    Operates by calling rotatePoint() on each x,y pair in turn. X and Y must
    have the same dimensions.

    Args:
        X (list): A list (or numpy array) of x coordinate values.
        Y (list): A list (or numpy array) of y coordinate values.
        rotation (float): The angle in degrees for the contour to be rotated by.
        centroid (tuple): A tuple containing the x,y coordinates of the centroid to rotate the contour about.

    Returns:
        tuple: A tuple containing a list of x coordinates and a list of y
        coordinates.
    '''

    rxs = []
    rys = []

    for nx, ny in zip(X, Y):
        rx, ry = rotatePoint((nx, ny), centroid, rotation)
        rxs.append(rx)
        rys.append(ry)

    return rxs, rys


def NormContour(X, Y, rawCentroid):
    '''
    Normalize the coordinates which make up a contour.

    Rescale the coordinates to values between 0 and 1 in both the x and y
    directions. The normalizing is performed using x or y width of the minimum
    bounding rectangle of the contour, whichever is largest. X and Y must have
    the same dimensions.

    Args:
        X (list): A list (or numpy array) of x coordinate values.
        Y (list): A list (or numpy array) of y coordinate values.
        rawCentroid (tuple): A tuple containing the x,y coordinates of the centroid of the contour.

    Returns:
        tuple: A tuple containing a list of normalized x coordinates, a list of
        normalized y coordinate and the normalized centroid.
    '''

    # find longest axis of rotated shape
    xwidth, ywidth, xmin, ymin = getBBoxDimensions(X, Y)
    if (xwidth > ywidth):
        normshape = xwidth
    elif (ywidth >= xwidth):
        normshape = ywidth

    norm_x = [(value - xmin) / normshape for value in X]
    norm_y = [(value - ymin) / normshape for value in Y]

    centroid = ((rawCentroid[0] - xmin) / normshape,
                (rawCentroid[1] - ymin) / normshape)

    return norm_x, norm_y, centroid


def CloseContour(X, Y):
    '''
    Close an opened polygon.

    Args:
        X (list): A list (or numpy array) of x coordinate values.
        Y (list): A list (or numpy array) of y coordinate values.

    Returns:
        tuple: A tuple containing the X and Y lists of coordinates where the
        first and last elements are equal.
    '''
    if ((X[0] != X[-1]) or (Y[0] != Y[-1])):
        X = X + [X[0]]
        Y = Y + [Y[0]]

    return X, Y


def ContourArea(X, Y):
    '''
    Compute the area of an irregular polygon.

    Ensures the contour is closed before processing, but does not modify
    X or Y outside the scope of this method. Algorithm taken from
    http://paulbourke.net/geometry/polygonmesh/.

    Args:
        X (list): A list (or numpy array) of x coordinate values.
        Y (list): A list (or numpy array) of y coordinate values.

    Returns:
        float: The area of the input polygon.
    '''

    # Check the contour provided is closed
    X, Y = CloseContour(X, Y)

    Sum = 0.

    for i in xrange(len(X) - 1):
        Sum += (X[i] * Y[i + 1]) - (X[i + 1] * Y[i])

    return abs(0.5 * Sum)


def ContourCentroid(X, Y):
    '''
    Compute the centroid of an irregular polygon.

    Ensures the contour is closed before processing, but does not modify
    X or Y outside the scope of this method. Algorithm taken from
    http://paulbourke.net/geometry/polygonmesh/.

    Args:
        X (list): A list (or numpy array) of x coordinate values.
        Y (list): A list (or numpy array) of y coordinate values.

    Returns:
        tuple: A tuple containing the (x,y) coordinate of the center of the
        input polygon.
    '''

    # Check the contour provided is closed
    X, Y = CloseContour(X, Y)

    Area = ContourArea(X, Y)

    Cx = 0.
    Cy = 0.

    for i in xrange(len(X) - 1):
        const = (X[i] * Y[i + 1]) - (X[i + 1] * Y[i])

        Cx += (X[i] + X[i + 1]) * const
        Cy += (Y[i] + Y[i + 1]) * const

    AreaFactor = (1. / (6. * Area))

    Cx *= AreaFactor
    Cy *= AreaFactor

    return (abs(Cx), abs(Cy))


def CalculateEFD(X, Y, harmonics=10):
    '''
    Compute the Elliptical Fourier Descriptors for a polygon.

    Implements Kuhl and Giardina method of computing the coefficients
    An, Bn, Cn, Dn for a specified number of harmonics. This code is adapted
    from the pyefd module. See the original paper for more detail:

    Kuhl, FP and Giardina, CR (1982). Elliptic Fourier features of a closed
    contour. Computer graphics and image processing, 18(3), 236-258.

    Args:
        X (list): A list (or numpy array) of x coordinate values.
        Y (list): A list (or numpy array) of y coordinate values.
        harmonics (int): The number of harmonics to compute for the given shape, defaults to 10.

    Returns:
        numpy.ndarray: A numpy array of shape (harmonics, 4) representing the
        four coefficients for each harmonic computed.
    '''

    contour = np.array([(x, y) for x, y in zip(X, Y)])

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
        a_n = const * np.sum((dxy[:, 1] / dt) * d_cos_phi_n)
        b_n = const * np.sum((dxy[:, 1] / dt) * d_sin_phi_n)
        c_n = const * np.sum((dxy[:, 0] / dt) * d_cos_phi_n)
        d_n = const * np.sum((dxy[:, 0] / dt) * d_sin_phi_n)
        coeffs[n - 1, :] = a_n, b_n, c_n, d_n

    return coeffs


def inverse_transform(coeffs, locus=(0., 0.), n=300, harmonic=10):
    '''
    Perform an inverse fourier transform to convert the coefficients back into
    spatial coordinates.

    Implements Kuhl and Giardina method of computing the performing the
    transform for a specified number of harmonics. This code is adapted
    from the pyefd module. See the original paper for more detail:

    Kuhl, FP and Giardina, CR (1982). Elliptic Fourier features of a closed
    contour. Computer graphics and image processing, 18(3), 236-258.

    Args:
        coeffs (numpy.ndarray): A numpy array of shape (n, 4) representing the four coefficients for each harmonic computed.
        locus (tuple): The x,y coordinates of the centroid of the contour being generated. Use calculate_dc_coefficients() to generate the correct locus for a shape.
        n (int): The number of coordinate pairs to compute. A larger value will result in a more complex shape at the expense of increased computational time. Defaults to 300.
        harmonics (int): The number of harmonics to be used to generate coordinates, defaults to 10. Must be <= coeffs.shape[0]. Supply a smaller value to produce coordinates for a more generalized shape.

    Returns:
        numpy.ndarray: A numpy array of shape (harmonics, 4) representing the
        four coefficients for each harmonic computed.
    '''

    t = np.linspace(0, 1., n)
    xt = np.ones((n,)) * locus[0]
    yt = np.ones((n,)) * locus[1]

    for n in xrange(harmonic):

        xt += ((coeffs[n, 2] * np.cos(2. * (n + 1.) * np.pi * t)) +
               (coeffs[n, 3] * np.sin(2. * (n + 1.) * np.pi * t)))

        yt += ((coeffs[n, 0] * np.cos(2. * (n + 1.) * np.pi * t)) +
               (coeffs[n, 1] * np.sin(2. * (n + 1.) * np.pi * t)))

        if n == harmonic - 1:
            return xt, yt


def InitPlot():
    '''
    Set up the axes for plotting, ensuring that x and y dimensions are equal.

    Returns:
        matplotlib.axes.Axes: Matplotlib axis instance.
    '''
    ax = plt.gca()
    ax.axis('equal')

    return ax


def PlotEllipse(ax, x, y, color='k', width=1.):
    '''
    Plots an ellipse represented as a series of x and y coordinates on a given
    axis.

    Args:
        ax (matplotlib.axes.Axes): Matplotlib axis instance.
        x (list): A list (or numpy array) of x coordinate values.
        y (list): A list (or numpy array) of y coordinate values.
        color (string): A matplotlib color string to color the line used to plot the ellipse. Defaults to k (black).
        width (float): The width of the plotted line. Defaults to 1.
    '''
    ax.plot(x, y, color, linewidth=width)


def SavePlot(ax, harmonic, filename, figformat='png'):
    '''
    Wrapper around the savefig method.

    Call this method to add a title identifying the harmonic being plotted, and
    save the plot to a file. Note that harmonic is simply an int value to be
    appended to the plot title, it does not select a harmonic to plot.

    The figformat argumet can take any value which matplotlib understands, which
    varies by system. To see a full list suitable for your matplotlib instance,
    call plt.gcf().canvas.get_supported_filetypes().

    Args:
        ax (matplotlib.axes.Axes): Matplotlib axis instance.
        harmonic (int): The harmonic which is being plotted.
        filename (string): A complete path and filename, without an extension, for the saved plot.
        figformat (string): A string denoting the format to save the figure as. Defaults to png.

    '''
    ax.set_title('Harmonic: {0}'.format(harmonic))
    plt.savefig('{0}_{1}.{2}'.format(filename, harmonic, figformat))
    plt.clf()


def PlotContour(ax, contour, color='b', width=1.):
    '''
    Plots a contour on a given axis.

    Args:
        ax (matplotlib.axes.Axes): Matplotlib axis instance.
        contour (numpy.ndarray): A numpy array of shape (n, 2) representing the input contour.
        color (string): A matplotlib color string to color the line used to plot the contour. Defaults to b (blue).
        width (float): The width of the plotted line. Defaults to 1.
    '''
    ax.plot(contour[:, 0], contour[:, 1], color, linewidth=width)


def AverageCoefficients(coeffList):
    '''
    Average the coefficients contained in the list of coefficient arrays,
    coeffList.

    This method is outlined in:

    2-D particle shape averaging and comparison using Fourier descriptors:
    Powder Technology Volume 104, Issue 2, 1 September 1999, Pages 180-189

    Args:
        coeffList (list): A list of coefficient arrays to be averaged.

    Returns:
        numpy.ndarray: A numpy array containing the average An, Bn, Cn, Dn
        coefficient values.
    '''

    nHarmonics = coeffList[0].shape[0]
    coeffsum = np.zeros((nHarmonics, 4))

    for coeff in coeffList:
        coeffsum += coeff

    coeffsum /= float(len(coeffList))

    return coeffsum


def AverageSD(coeffList, avgcoeffs):
    '''
    Use the coefficients contained in the list of coefficient arrays,
    coeffList, and the average coefficient values to compute the standard
    deviation of series of ellipses.

    This method is outlined in:

    2-D particle shape averaging and comparison using Fourier descriptors:
    Powder Technology Volume 104, Issue 2, 1 September 1999, Pages 180-189

    Args:
        coeffList (list): A list of coefficient arrays to be averaged.
        avgcoeffs (numpy.ndarray): A numpy array containing the average coefficient values, generated by calling AverageCoefficients().

    Returns:
        numpy.ndarray: A numpy array containing the standard deviation
        An, Bn, Cn, Dn coefficient values.
    '''
    nHarmonics = avgcoeffs.shape[0]
    coeffsum = np.zeros((nHarmonics, 4))

    for coeff in coeffList:
        coeffsum += (coeff ** 2.)

    return (coeffsum / float(len(coeffList) - 1)) - (avgcoeffs ** 2.)


def Nyquist(contour):
    '''
    Returns the maximum number of harmonics that can be computed for a given
    contour, the nyquist freqency.

    See this paper for details:
    C. Costa et al. / Postharvest Biology and Technology 54 (2009) 38-47

    Args:
        contour (numpy.ndarray): A numpy array of shape (n, 2) representing the input contour.

    Returns:
        int: The nyquist frequency, expressed as a number of harmonics.
    '''
    return contour.shape[0] / 2


def FourierPower(coeffs, contour, threshold=0.9999):
    '''
    Compute the total Fourier power and find the minium number of harmonics
    required to exceed the threshold fraction of the total power.

    This is a good method for identifying the number of harmonics to use to
    describe a polygon. For more details see:

    C. Costa et al. / Postharvest Biology and Technology 54 (2009) 38-47

    Warning:
        The number of coeffs must be >= the nyquist freqency.

    Args:
        coeffs (numpy.ndarray): A numpy array of shape (n, 4) representing the four coefficients for each harmonic computed.
        contour (numpy.ndarray): A numpy array of shape (n, 2) representing the input contour.
        threshold (float): The threshold fraction of the total Fourier power, the default is 0.9999.

    Returns:
        int: The number of harmonics required to represent the contour above the
        threshold Fourier power.

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

        if (currentPower / totalPower) > threshold:
            return i + 1


def normalize_efd(coeffs, size_invariant=True):
    '''
    Normalize the Elliptical Fourier Descriptor coefficients for a polygon.

    Implements Kuhl and Giardina method of normalizing the coefficients
    An, Bn, Cn, Dn. Performs 3 separate normalizations. First, it makes the data
    location invariant by re-scaling the data to a common origin. Secondly, the
    data is rotated with respect to the major axis. Thirdly, the coefficients
    are normalized with regard to the absolute value of A_1. This code is
    adapted from the pyefd module. See the original paper for more detail:

    Kuhl, FP and Giardina, CR (1982). Elliptic Fourier features of a closed
    contour. Computer graphics and image processing, 18(3), 236-258.

    Args:
        coeffs (numpy.ndarray): A numpy array of shape (n, 4) representing the four coefficients for each harmonic computed.
        size_invariant (bool): Set to True (the default) to perform the third normalization and false to return the data withot this processing step.

    Returns:
        numpy.ndarray: A numpy array of shape (harmonics, 4) representing the
        four coefficients for each harmonic computed.
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


def calculate_dc_coefficients(X, Y):
    '''
    Compute the dc coefficients, used as the locus when calling
    inverse_transform().

    This code is adapted from the pyefd module. See the original paper for
    more detail:

    Kuhl, FP and Giardina, CR (1982). Elliptic Fourier features of a closed
    contour. Computer graphics and image processing, 18(3), 236-258.

    Args:
        X (list): A list (or numpy array) of x coordinate values.
        Y (list): A list (or numpy array) of y coordinate values.

    Returns:
        tuple: A tuple containing the c and d coefficients.

    Todo: Test this against other efd modules to check that the x,y to y,x
    conversion is correct.
    '''

    contour = np.array([(x, y) for x, y in zip(X, Y)])

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
    return (contour[0, 0] + A0, contour[0, 1] + C0)


def LoadGeometries(filename):
    '''
    Takes a filename and uses pyshp to load it, returning a list of
    shapefile.ShapeRecord instances.

    This list can be iterated over, passing the individual shape instances
    to ProcessGeometry() one by one. There is no input handling if a non-polygon
    shapefile is passed in, that will result in undefined behavior.

    Args:
        filename (string): A filename with optional full path pointing to an  ESRI shapefile to be loaded by the pyshp module. The file extension is optional.

    Returns:
        list: A list of shapefile._ShapeRecord objects representing each polygon
        geometry in the shapefile.
    '''
    shp = sf.Reader(filename)
    return shp.shapeRecords()


def ProcessGeometry(shape):
    '''
    Method to handle all the geometry conversion that may be needed by the rest
    of the EFD code.

    Method which takes a single shape instance from a shapefile
    eg shp.Reader('shapefile.shp').shapeRecords()[n]
    where n is the index of the shape within a multipart geometry. This results
    in the contour, coordinate list and centroid data computed for the input
    polygon being normalized and returned to the user.

    Args:
        shapefile._ShapeRecord: A shapefile object representing the geometry and attributes of a single polygon from a multipart shapefile.

    Returns:
        tuple: A tuple containing a list of normalized x coordinates, a list of
        normalized y coordinates, contour (a list of [x,y] coordinate pairs,
        normalized about the shape's centroid) and the normalized coordinate
        centroid.
    '''
    x = []
    y = []

    for point in shape.shape.points:
        x.append(point[0])
        y.append(point[1])

    centroid = ContourCentroid(x, y)
    X, Y, NormCentroid = NormContour(x, y, centroid)

    contour = np.array([(a, b) for a, b in zip(X, Y)])

    return X, Y, contour, NormCentroid


def rotatePoint(point, centerPoint, angle):
    '''
    Rotates a point counter-clockwise around centerPoint.

    The angle to rotate by is supplied in degrees. Code based on:
    https://gist.github.com/somada141/d81a05f172bb2df26a2c

    Args:
        point (tuple): The point to be rotated, represented as an (x,y) tuple.
        centerPoint (tuple): The point to be rotated about, represented as an(x,y) tuple.
        angle (float): The angle to rotate point by, in the counter-clockwise direction.

    Returns:
        tuple: A tuple representing the rotated point, (x,y).
    '''
    angle = np.radians(angle)
    temp_point = point[0] - centerPoint[0], point[1] - centerPoint[1]
    temp_point = (temp_point[0] * np.cos(angle) - temp_point[1] *
                  np.sin(angle), temp_point[0] * np.sin(angle) +
                  temp_point[1] * np.cos(angle))

    temp_point = temp_point[0] + centerPoint[0], temp_point[1] + centerPoint[1]
    return temp_point[0], temp_point[1]


def getBBoxDimensions(x, y):
    '''
    Returns the width in the x and y dimensions and the maximum x and y
    coordinates for the bounding box of a given list of x and y coordinates.

    Args:
        x (list): A list (or numpy array) of x coordinate values.
        y (list): A list (or numpy array) of y coordinate values.
    Returns:
        tuple: A four-tuple representing (width in the x direction, width in the
        y direction, the minimum x coordinate and the minimum y coordinate).
    '''
    xmin = float(min(x))
    ymin = float(min(y))

    xmax = float(max(x))
    ymax = float(max(y))

    return xmax - xmin, ymax - ymin, xmin, ymin
