# -*- coding: utf-8 -*-
import glob
import pickle
from sklearn.externals import joblib
from unittest import TestCase


class LibRegressorsTestCase(TestCase):

    def test_all_c1_n1_or_n1_n9(self):
        regressor = joblib.load("restraintlib/lib/regressors/All-C1'-N1 or C1'-N9.joblib")
        value, sigma = regressor.predict([[10.0]], return_std=True)
        self.assertAlmostEqual(value[0], 1.4874497)
        self.assertAlmostEqual(sigma[0], 0.011119879)

        with open("restraintlib/lib/regressors/All-C1'-N1 or C1'-N9.pickle", 'r') as p_file:
            regressor = pickle.load(p_file)
        value, sigma = regressor.predict([[10.0]], return_std=True)
        self.assertAlmostEqual(value[0], 1.4874497)
        self.assertAlmostEqual(sigma[0], 0.011119879)

    def test_joblib_vs_pickle(self):
        names = glob.glob("restraintlib/lib/regressors/*.joblib")
        for name in names:
            print(name)
            regressor_joblib = joblib.load(name)
            with open(name.replace('joblib', 'pickle'), 'r') as p_file:
                regressor_pickle = pickle.load(p_file)

            for x in [5., 10., 22.5, 30., 60]:
                value_joblib, sigma_joblib = regressor_joblib.predict([[x]], return_std=True)
                value_pickle, sigma_pickle = regressor_pickle.predict([[x]], return_std=True)
                self.assertAlmostEqual(value_joblib[0], value_pickle[0])
                self.assertAlmostEqual(sigma_joblib[0], sigma_pickle[0])
