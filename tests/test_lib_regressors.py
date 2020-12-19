# -*- coding: utf-8 -*-
from unittest import TestCase

from restraintlib.restraints import load_function


class LibRegressorsTestCase(TestCase):

    def assertFunctional(self, function_name, x, expected_value, expected_sigma):
        regressor = load_function(function_name)
        value, sigma = regressor.predict([[x]], return_std=True)
        print('self.assertFunctional("%s", %.1f, %.7f, %.10f)' % (function_name, x, value[0], sigma[0]))
        self.assertAlmostEqual(value[0], expected_value)
        self.assertAlmostEqual(sigma[0], expected_sigma)

    def test_pickle(self):
        self.assertFunctional("All-C1'-N1 or C1'-N9.pickle", 10.0, 1.4874497, 0.011119879)
        self.assertFunctional("All-C1'-O4'.pickle", 10.0, 1.4095869, 0.0092497524)
        self.assertFunctional("deoxyribose-C2'-endo-C1'-C2'-C3'.pickle", 10.0, 108.46958654, 0.705017644)
        self.assertFunctional("deoxyribose-C2'-endo-C1'-O4'-C4'.pickle", 10.0, 113.33767771, 0.654150846)
        self.assertFunctional("deoxyribose-C2'-endo-C2'-C3'-C4'.pickle", 10.0, 107.58851579, 0.853301519)
        self.assertFunctional("deoxyribose-C2'-endo-C3'-C4'-O4'.pickle", 10.0, 107.49313064, 0.867046559)
        self.assertFunctional("deoxyribose-C3'-endo-C1'-C2'-C3'.pickle", 10.0, 109.23909162, 0.981169342)
        self.assertFunctional("deoxyribose-C3'-endo-C1'-O4'-C4'.pickle", 10.0, 112.37972166, 0.946601366)
        self.assertFunctional("deoxyribose-C3'-endo-C2'-C3'-C4'.pickle", 10.0, 106.43252974, 0.842729834)
        self.assertFunctional("deoxyribose-C3'-endo-C3'-C4'-O4'.pickle", 10.0, 109.46664887, 0.696293795)
        self.assertFunctional("deoxyribose-Other-C1'-C2'-C3'.pickle", 10.0, 107.0025446, 1.4534300232)
        self.assertFunctional("deoxyribose-Other-C1'-O4'-C4'.pickle", 10.0, 112.9288726, 1.3327459941)
        self.assertFunctional("deoxyribose-Other-C2'-C3'-C4'.pickle", 10.0, 106.9423613, 1.2974987305)
        self.assertFunctional("deoxyribose-Other-C3'-C4'-O4'.pickle", 10.0, 109.3494080, 1.7140458693)
        self.assertFunctional("purine-C1'-N1-C2 or C1'-N9-C4.pickle", 10.0, 126.9829216, 1.5048922785)
        self.assertFunctional("purine-C1'-N1-C6 or C1'-N9-C8.pickle", 10.0, 127.1839575, 1.6464197451)
        self.assertFunctional("purine-N1-C1'-C2' or N9-C1'-C2'.pickle", 10.0, 113.6571028, 1.2607817385)
        self.assertFunctional("purine-N1-C1'-O4' or N9-C1'-O4'.pickle", 10.0, 107.3662676, 1.1939858801)
        self.assertFunctional("pyrimidine-C1'-N1-C2 or C1'-N9-C4.pickle", 10.0, 117.9132730, 1.3865739592)
        self.assertFunctional("pyrimidine-C1'-N1-C6 or C1'-N9-C8.pickle", 10.0, 120.7343454, 1.2841006880)
        self.assertFunctional("pyrimidine-N1-C1'-C2' or N9-C1'-C2'.pickle", 10.0, 113.5856496, 1.3802839915)
        self.assertFunctional("pyrimidine-N1-C1'-O4' or N9-C1'-O4'.pickle", 10.0, 108.0619242, 1.1917431290)
        self.assertFunctional("ribose-C2'-endo-C1'-C2'-C3'.pickle", 10.0, 107.6786200, 0.8340928962)
        self.assertFunctional("ribose-C2'-endo-C1'-O4'-C4'.pickle", 10.0, 113.6515472, 0.8514767462)
        self.assertFunctional("ribose-C2'-endo-C2'-C3'-C4'.pickle", 10.0, 106.5254955, 0.8538892370)
        self.assertFunctional("ribose-C2'-endo-C3'-C4'-O4'.pickle", 10.0, 109.0108704, 0.8917961243)
        self.assertFunctional("ribose-C3'-endo-C1'-C2'-C3'.pickle", 10.0, 106.5080838, 1.1667549808)
        self.assertFunctional("ribose-C3'-endo-C1'-O4'-C4'.pickle", 10.0, 114.2986980, 0.7747160588)
        self.assertFunctional("ribose-C3'-endo-C2'-C3'-C4'.pickle", 10.0, 108.9975249, 0.9434856130)
        self.assertFunctional("ribose-C3'-endo-C3'-C4'-O4'.pickle", 10.0, 109.3762512, 0.8818958553)
        self.assertFunctional("ribose-Other-C1'-C2'-C3'.pickle", 10.0, 108.1068697, 2.3693564847)
        self.assertFunctional("ribose-Other-C1'-O4'-C4'.pickle", 10.0, 116.0915706, 1.6088767710)
        self.assertFunctional("ribose-Other-C2'-C3'-C4'.pickle", 10.0, 106.1877587, 1.4164389140)
        self.assertFunctional("ribose-Other-C3'-C4'-O4'.pickle", 10.0, 104.7131548, 1.0766946431)
