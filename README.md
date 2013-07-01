Some experiments for deriving the location of a cell tower based signal strength
readings from devices at known locations.

The methods implemented here are very similar to the one presented in
[this StackExchange thread][link-stackexchange], which consists of a nonlinear
least-squares minimization where data points are weighted on the inverse of
their distance to the cell tower.

I have also included an implementation of the [RANSAC algorithm][ransac-wikipedia]
imported from the repository at [FredrikAppelros/ransac][ransac-repository].
This algorithm can be used to filter outliers from the input data.

The directory data/ contains sample measures used to compare different
variations of the method implemented in run.py. To run them, use:

    python run.py data/measures-262-1-34823-45550.csv

You should see a table comparing results from using the original method with and
without weighting, and with and without running RANSAC before the actual
optimization.

padnums.py is an adapted version of the found module at [this blog post][link-padnums]

Dependencies: scipy and numpy.

[ransac-repository]: https://github.com/FredrikAppelros/ransac
[ransac-wikipedia]: https://en.wikipedia.org/wiki/RANSAC
[link-stackexchange]: http://gis.stackexchange.com/questions/40660/trilateration-algorithm-for-n-amount-of-points?lq=1
[link-padnums]: http://ginstrom.com/scribbles/2007/09/04/pretty-printing-a-table-in-python/
