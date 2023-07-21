import multiprocessing as mp
import os
import time
from collections import Counter
from os import path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from bbmix.models import MixtureBinomialSparseBatch


###old model, do not use now


class Metaquad():
    def __init__(self, AD, DP, variant_names=None, dataset_name=None):

        self.ad = AD
        self.dp = DP

        if AD.shape[0] != DP.shape[0]:
            print('AD and DP length do not match!')
        else:
            print(str(AD.shape[0]) + ' variants detected...')

        if variant_names is not None:
            #  check length
            if len(variant_names) != self.ad.shape[0]:
                print('No. of variant names does not match length of AD!')
            else:
                self.variants = variant_names
                print("Variant names detected...")

        else:
            self.variants = None

        if dataset_name is not None:
            self.dataset = dataset_name

    def _check_outdir_exist(self, out_dir):
        if path.exists(out_dir) is not True:
            try:
                os.mkdir(out_dir)
                return True
            except:
                print("Can't make directory, do you have permission?")
                return False
        else:
            print('Out directory already exists, overwriting content inside...')
            return True

    def _addVariantNames(self, valid_rows):
        var_col = 'variant_name'
        variants = np.array(self.variants)
        df = self.df
        df_zeros = pd.DataFrame(np.zeros((variants.shape[0] - df.shape[0], df.shape[1])), columns=df.columns)
        df[var_col] = variants[valid_rows]
        df_zeros[var_col] = list(set(variants) - set(df[var_col]))
        self.df = pd.concat([df, df_zeros], axis=0, ignore_index=True)

    def _batchify(self, batch_size, x, y, valid_row_sizes):
        n_samples = valid_row_sizes.shape[0]
        start, item_start = 0, 0
        while start < n_samples:
            end = min(batch_size + start, n_samples)
            _row_sizes = valid_row_sizes[start:end]

            item_end = item_start + np.sum(_row_sizes)
            _x, _y = x[item_start:item_end], y[item_start:item_end]
            _ad = np.squeeze(np.asarray(self.ad[_x, _y]))
            _dp = np.squeeze(np.asarray(self.dp[_x, _y]))

            start = end
            item_start = item_end
            yield (_ad, _dp, _row_sizes)
        assert (item_start == x.shape[0])


    def fit_deltaBIC(self, out_dir, minDP=10, minAD=1, export_csv=True, nproc=30, batch_size=128):

        n_variants = self.dp.shape[0]

        #Batch size adjustment
        adj_bs = min(n_variants // nproc, batch_size)
        print('CPUs used: {}, batch size: {} {}'.format(nproc,
                                                        adj_bs,
                                                        "" if adj_bs == batch_size else "(adjusted to avoid idle processes)"))
        print("Fitting in sparse mode...")
        t0 = time.time()

        print("Initializing fit(mode: deltaBIC) on " + str(self.ad.shape[0]) + " variants...")

        dp_row, dp_col = np.nonzero(self.dp >= minDP)
        ad_row, _ = np.nonzero(self.ad >= minAD)

        # only variant with at least one valid dp and ad records are included
        valid_rows = np.intersect1d(dp_row, ad_row)
        # filter variants
        x, y = zip(*[(r, c) for r, c in zip(dp_row, dp_col) if r in valid_rows])

        # split batch
        x, y = np.array(x), np.array(y)
        valid_row_sizes = Counter(x)
        valid_row_sizes = np.array([valid_row_sizes[r_idx] for r_idx in valid_rows])

        assert (np.sum(valid_row_sizes) == x.shape[0])

        with mp.Pool(processes=nproc) as pool:
            results = pool.starmap_async(fit_batch, self._batchify(batch_size, x, y, valid_row_sizes)).get()
        self.df = pd.concat([pd.DataFrame(res) for res in results], axis=0, ignore_index=False)

        t1 = time.time()
        print("deltaBIC was calculated for " + str(self.ad.shape[0]) + " variants and took:%.2f minutes" % (
                    (t1 - t0) / 60))


        if self.variants is not None:
            self._addVariantNames(valid_rows)

        # sort csv
        self.sorted_df = self.df.sort_values(by=['deltaBIC'], ascending=False)

        if export_csv is True:
            if self._check_outdir_exist(out_dir) is True:
                self.sorted_df.to_csv(out_dir + '/BIC_params.csv', index=False)
            else:
                self.sorted_df.to_csv('BIC_params.csv', index=False)

        self.df.to_csv(out_dir + '/debug_unsorted_BIC_params.csv', index=False)
        # return csv
        return self.df

    def selectInformativeVariants(self, min_samples=5, export_heatmap=False, out_dir=None,cutoff=10):

        if self.df is None:
            print('No Fitted model!')
        else:
            if out_dir is not None:
                if path.exists(out_dir) is not True:
                    try:
                        os.mkdir(out_dir)
                    except:
                        print("ERROR")
            else:
                print('Output directory already exists, overwriting...')


            print('deltaBIC cutoff = ', cutoff)
            self.sorted_df['PASS_cutoff'] = self.sorted_df.deltaBIC.apply(lambda x: True if x >= cutoff else False)
            self.sorted_df['PASS_MINSAMPLES'] = self.sorted_df.num_samples_minor_cpt.apply(
                lambda x: True if x >= min_samples else False)

            self.final_df = self.sorted_df[(self.sorted_df.PASS_cutoff == True) & (self.sorted_df.PASS_MINSAMPLES == True)]


            print('Number of variants passing threshold: ' + str(len(self.final_df['variant_name'])))

            if len(self.final_df['variant_name']) != 0:
                passed_variants = self.final_df['variant_name']
                idx = [self.variants.index(i) for i in passed_variants]

                best_ad = self.ad[idx]
                best_dp = self.dp[idx]
            else:
                print(
                    "No informative variants detected!")

        self.sorted_df.to_csv(out_dir + '/BIC_params.csv', index=False)


        plt.plot(self.sorted_df.index.sort_values(ascending=False),self.sorted_df.deltaBIC)
        plt.axhline(y=cutoff, color="black", linestyle='--', label="cutoff")
        plt.legend()
        plt.ylabel("\u0394BIC")
        plt.xlabel("Cumulative number of mutations")
        plt.savefig(out_dir + '/deltaBIC_cdf.pdf')




        if self.variants is not None:
            renamed_vars = []
            for var in passed_variants:
                renamed_vars.append((var.split('_')[1] + var.split('_')[2] + '>' + var.split('_')[3]))

            with open(out_dir + '/passed_variant_names.txt', "w+") as var_file:
                var_file.write('\n'.join(str(var) for var in renamed_vars))

        if export_heatmap is True:
            af = best_ad / best_dp
            fig, ax = plt.subplots(figsize=(15, 10))
            plt.title("Allele frequency of top variants")
            plt.style.use('seaborn-dark')
            if self.variants is not None:
                sns.heatmap(af, cmap='Greens', yticklabels=renamed_vars)
            else:
                sns.heatmap(af, cmap='Greens')
            plt.savefig(out_dir + '/top variants heatmap.pdf')

        return best_ad, best_dp

    def readParams(self, file):
        self.df = pd.read_csv(file)
        self.sorted_df = self.df.sort_values(by=['deltaBIC'], ascending=False)

        return self.df, self.sorted_df


def sparseMixBinFit(valid_ad, valid_dp, valid_row_sizes, fix_seed=None):
    # delta BIC
    if fix_seed is not None:
        np.random.seed(fix_seed)

    model1 = MixtureBinomialSparseBatch(n_components=1, tor=1e-20)
    params1 = model1.fit((valid_ad, valid_dp), valid_row_sizes, max_iters=500, early_stop=True)

    model2 = MixtureBinomialSparseBatch(n_components=2, tor=1e-20)
    params2 = model2.fit((valid_ad, valid_dp), valid_row_sizes, max_iters=500, early_stop=True)

    model3 = MixtureBinomialSparseBatch(n_components=5, tor=1e-20)
    params3 = model3.fit((valid_ad, valid_dp), valid_row_sizes, max_iters=500, early_stop=True)

    delta_BIC = model1.model_scores["BIC"] - model2.model_scores["BIC"]
    delta_BIC2= model1.model_scores["BIC"] - model3.model_scores["BIC"]

    p = params2[:, [0, 1]]
    pi = params2[:, [2, 3]]
    fraction_b_allele = np.min(p, axis=1) * np.array([pi[ith, idx] for ith, idx in enumerate(np.argmin(p, axis=1))])




    minor_cpt_n = np.min(pi, axis=1) * valid_row_sizes

    results = {
        "num_samples": valid_row_sizes,
        'deltaBIC2': delta_BIC2,
        'deltaBIC': delta_BIC,
        'params1': params1.tolist(),
        'params2': params2.tolist(),
        'model1BIC': model1.model_scores["BIC"],
        'model2BIC': model2.model_scores["BIC"],
        'fraction_b_allele': fraction_b_allele,
        'num_samples_minor_cpt': minor_cpt_n,
    }
    return results


def fit_batch(valid_ad, valid_dp, valid_row_sizes):
    basic_stats = basicStats(valid_ad, valid_dp, valid_row_sizes)
    results = sparseMixBinFit(valid_ad, valid_dp, valid_row_sizes)
    results.update(basic_stats)
    return results


def basicStats(valid_ad, valid_dp, valid_row_sizes):

    left, batch_size = 0, len(valid_row_sizes)
    stats = ['total_DP', 'median_DP', 'total_AD', 'median_AD', 'num_samples_nonzero_AD']
    batch_res = {name: np.empty(batch_size) for name in stats}

    for ith, smp_sz in enumerate(valid_row_sizes):
        right = left + smp_sz
        _d = valid_dp[left:right]
        _a = valid_ad[left:right]
        batch_res[stats[0]][ith] = np.sum(_d)
        batch_res[stats[1]][ith] = np.median(_d)
        batch_res[stats[2]][ith] = np.sum(_a)
        batch_res[stats[3]][ith] = np.median(_a)
        batch_res[stats[4]][ith] = np.count_nonzero(_a)

        left = right

    return batch_res

    total_DP = np.sum(_d, axis=1)
    median_DP = np.median(_d, axis=1)
    total_AD = np.sum(_a, axis=1)
    median_AD = np.median(_a, axis=1)
    non_zero = np.count_nonzero(_a, axis=1)
    return {'total_DP': total_DP,
            'median_DP': median_DP,
            'total_AD': total_AD,
            'median_AD': median_AD,
            'num_samples_nonzero_AD': non_zero
            }
