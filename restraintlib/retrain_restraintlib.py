import math
import os
import pickle

import pandas as pd
from sklearn.linear_model import BayesianRidge
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import WhiteKernel, ExpSineSquared
from sklearn.metrics import mean_squared_error, median_absolute_error, r2_score, mean_absolute_error

__author__ = "Marcin Kowiel, Dariusz Brzezinski"


def serialize_pickle(model, filename):
    # if python3 use protocol 4 or highest, to be compatible with python2 we can use maximally protocol 2
    print('Saving to', filename)
    with open(filename, 'wb') as pickle_file:
        pickle.dump(model, pickle_file, protocol=2)


def save_regressor_stats(regressor, x, y, group, y_col, stats_df):
    y_hat = regressor.predict(x)

    try:
        decision_func = "{0} + {1:.3f}".format(str(regressor.kernel_), regressor._y_train_mean)
    except:
        decision_func = "{0:.3f}x + {1:.3f}".format(regressor.coef_[0], regressor.intercept_)

    stats = pd.DataFrame({
        "Subgroup": [group],
        "Measurement": [y_col],
        "Coefficients": [decision_func],
        "R^2": [r2_score(y, y_hat)],
        "RMSE": [math.sqrt(mean_squared_error(y, y_hat))],
        "MAE": [mean_absolute_error(y, y_hat)],
        "MAD": [median_absolute_error(y, y_hat)]
    })

    return stats_df.append(stats)


def create_linear_regressors(df, x_col, y_cols, stats_df, dir_name='lib/regressors/'):
    if not os.path.exists(dir_name):
        os.makedirs(dir_name)

    for y_col in y_cols:
        for sugar in ["ribose", "deoxyribose"]:
            for conformation in ["C2'-endo", "C3'-endo", "Other"]:
                name = sugar + "-" + conformation + "-" + y_col.replace("/", " or ")
                x = df.loc[:, x_col][df.Sugar == sugar][df.Conformation == conformation]
                x = x.values.reshape(-1, 1)
                y = df.loc[:, y_col][df.Sugar == sugar][df.Conformation == conformation]

                gpr = BayesianRidge(normalize=True)
                gpr.fit(x, y)

                stats_df = save_regressor_stats(gpr, x, y, sugar + "-" + conformation, y_col, stats_df)
                serialize_pickle(gpr, os.path.join(dir_name, name + ".pickle"))

    return stats_df


def create_sine_regressors(df, x_col, y_cols, period, stats_df, use_base=True, dir_name='lib/regressors/'):
    if not os.path.exists(dir_name):
        os.makedirs(dir_name)

    for y_col in y_cols:
        if use_base:
            for base in ["purine", "pyrimidine"]:
                name = base + "-" + y_col.replace("/", " or ")
                x = df.loc[:, x_col][df.Base == base]
                x = x.values.reshape(-1, 1)
                y = df.loc[:, y_col][df.Base == base]

                gpr = GaussianProcessRegressor(
                    kernel=(ExpSineSquared(periodicity_bounds=(period, period)) +
                            WhiteKernel(noise_level_bounds=(1e-7, 1e7))),
                    n_restarts_optimizer=100,
                    random_state=23,
                    normalize_y=True
                )
                gpr.fit(x, y)

                stats_df = save_regressor_stats(gpr, x, y, base, y_col, stats_df)
                serialize_pickle(gpr, os.path.join(dir_name, name + ".pickle"))
        else:
            name = 'All-' + y_col.replace("/", " or ")
            x = df.loc[:, x_col]
            x = x.values.reshape(-1, 1)
            y = df.loc[:, y_col]

            gpr = GaussianProcessRegressor(
                kernel=(ExpSineSquared(periodicity_bounds=(period, period)) +
                        WhiteKernel(noise_level_bounds=(1e-7, 1e7))),
                n_restarts_optimizer=100,
                random_state=23,
                normalize_y=True
            )
            gpr.fit(x, y)

            stats_df = save_regressor_stats(gpr, x, y, "All", y_col, stats_df)
            serialize_pickle(gpr, os.path.join(dir_name, name + ".pickle"))

    return stats_df


def run():
    abs_data_path = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'data', 'combined_results.csv')
    abs_pickle_dir = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'lib', 'regressors')
    sugar_measurements = pd.read_csv(abs_data_path)
    print(sugar_measurements.head())
    stats_df = pd.DataFrame()
    stats_df = create_linear_regressors(sugar_measurements, "T_max",
                                        ["C1'-C2'-C3'", "C2'-C3'-C4'", "C3'-C4'-O4'", "C1'-O4'-C4'"], 
                                        stats_df,
                                        dir_name=abs_pickle_dir)
    stats_df = create_sine_regressors(sugar_measurements, "TCHI", ["C1'-N1/C1'-N9", "C1'-O4'"], 180, stats_df,
                                      use_base=False, dir_name=abs_pickle_dir)
    stats_df = create_sine_regressors(sugar_measurements, "TCHI",
                                      ["N1-C1'-C2'/N9-C1'-C2'", "C1'-N1-C2/C1'-N9-C4", "C1'-N1-C6/C1'-N9-C8",
                                       "N1-C1'-O4'/N9-C1'-O4'"], 360, stats_df, dir_name=abs_pickle_dir)
    print(stats_df.head())


if __name__ == "__main__":
    run()
