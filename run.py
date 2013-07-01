import sys
import csv

import padnums
from ransac import ransac

import numpy as np
import scipy.optimize as optimize

class NonLinearLeastSquaresModel(ransac.Model):
    def __init__(self, weighted = False, initial_estimate = None):
        self.weighted = weighted
        self.initial_estimate = initial_estimate

        if self.weighted:
            self.residuals = self.residuals_weighted
        else:
            self.residuals = self.residuals_unweighted

    def weight(self, point):
        return 1.0 / point[-1]**2

    def residuals_weighted(self, current):
        return np.array(
                [self.weight(d) * np.linalg.norm(d[:-1] - current)
                for d in self.data])

    def residuals_unweighted(self, current):
        return np.array([np.linalg.norm(d[:-1] - current) for d in self.data])

    def fit(self, data):
        self.data = data

        if self.initial_estimate is None:
            self.initial_estimate = data[0][:-1]

        (x, y), status = optimize.leastsq(self.residuals, self.initial_estimate)
        if status not in (1, 2, 3, 4):
            raise ValueError("Failed to find least squares solution!")

        self.params = (x, y)
        self.residual = sum(r**2 for r in self.residuals(np.array(self.params)))

    def distance(self, data):
        self.data = data
        return self.residuals(self.params)

def load_data(csvpath):
    data = []
    with open(csvpath, "r") as cvsf:
        csvrd = csv.reader(cvsf, delimiter = ";")

        for row in csvrd:
            latitude, longitude, asu = map(float, row)
            if asu > 0: # filter out really bad readings
                data.append((latitude, longitude, asu))

    return data

if __name__ == "__main__":
    data = load_data(sys.argv[1])
    actual_location = (49.754592, 9.961313)

    def get_error(result):
        return np.linalg.norm(np.array(result) - actual_location)

    table = [("Method", "Latitude", "Longitude", "Residual", "Error")]
    table.append(("Actual location", 49.754592, 9.961313, "--", "--"))

    model = NonLinearLeastSquaresModel()
    model.fit(data)

    table.append(("Nonlinear least squares (unweighted)",
        model.params[0], model.params[1], model.residual,
        get_error(model.params)))

    model = NonLinearLeastSquaresModel(weighted = True)
    model.fit(data)

    table.append(("Nonlinear least squares (weighted)",
        model.params[0], model.params[1], model.residual,
        get_error(model.params)))

    (latitude, longitude), inliers, error = ransac.ransac(
        data, NonLinearLeastSquaresModel(),
        iterations = 10, min_samples = 3, min_inliers = 0.6, eps = 1e-3)

    table.append(("RANSAC (unweighted, %d/%d points)" % (len(inliers), len(data)),
        latitude, longitude, model.residual, get_error(model.params)))

    (latitude, longitude), inliers, error = ransac.ransac(
        data, NonLinearLeastSquaresModel(weighted = True),
        iterations = 10, min_samples = 3, min_inliers = 0.6, eps = 1e-3)

    table.append(("RANSAC (weighted, %d/%d points)" % (len(inliers), len(data)),
        latitude, longitude, model.residual, get_error(model.params)))

    padnums.pprint_table(sys.stdout, table)
