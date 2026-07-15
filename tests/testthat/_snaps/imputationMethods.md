# min2impute

    Code
      D_min
    Output
              sample1   sample2   sample3
      pep_1  8.251585  7.719108  9.232288
      pep_2  7.918788  9.423502  8.084255
      pep_3  8.751157  7.716832  9.693175
      pep_4  7.550698 10.295697  9.282804
      pep_5  8.611378  7.604930  9.477185
      pep_6  9.508241 11.248634 11.012259
      pep_7  7.856285  8.314790  7.529904
      pep_8  7.512251  9.095206  8.592066
      pep_9  8.185470 10.489451 10.438342
      pep_10 8.459545  7.729362  9.325311

# colMeanImputation

    Code
      D_min
    Output
                   sample1.pep_1 sample1.pep_2 sample1.pep_3 sample1.pep_4
      sample1_run1      17.57193      15.83758      23.73870      18.29231
      sample1_run2      16.50317      24.09699      21.87331      15.10140
      sample1_run3      18.70053      19.72942      17.50231      21.22343
                   sample1.pep_5 sample1.pep_6 sample1.pep_7 sample1.pep_8
      sample1_run1      17.22276      19.01648      15.72499      15.02450
      sample1_run2      23.05993      20.58218      22.14794      21.34943
      sample1_run3      22.63031      19.72942      15.71257      19.18439
                   sample1.pep_9 sample1.pep_10 sample2.pep_1 sample2.pep_2
      sample1_run1      16.37094       16.91909      15.43822      23.27206
      sample1_run2      21.18414       19.92332      19.70244      18.84700
      sample1_run3      18.79415       24.08768      17.18089      24.87634
                   sample2.pep_3 sample2.pep_4 sample2.pep_5 sample2.pep_6
      sample1_run1      15.43366      24.89333      15.89846      23.17690
      sample1_run2      16.56922      20.59139      16.46224      23.13152
      sample1_run3      18.48462      19.43643      15.20986      22.49727
                   sample2.pep_7 sample2.pep_8 sample2.pep_9 sample2.pep_10
      sample1_run1      24.38314      21.20169      20.50695       20.86509
      sample1_run2      21.84061      18.87477      19.05310       15.45872
      sample1_run3      16.62958      18.19041      20.97890       20.87997
                   sample3.pep_1 sample3.pep_2 sample3.pep_3 sample3.pep_4
      sample1_run1      23.24674      16.16851      21.03924      18.56561
      sample1_run2      18.46458      20.08279      22.56324      19.81107
      sample1_run3      18.83131      18.90382      19.38635      19.63025
                   sample3.pep_5 sample3.pep_6 sample3.pep_7 sample3.pep_8
      sample1_run1      24.65195      24.67561      21.03924      18.32156
      sample1_run2      18.95437      22.02452      18.33825      17.18413
      sample1_run3      20.83005      18.83131      15.05981      17.18942
                   sample3.pep_9 sample3.pep_10
      sample1_run1      21.53844       21.14550
      sample1_run2      20.87668       19.81107
      sample1_run3      21.00017       18.65062

# missForest

    Code
      D_missForest
    Output
             sample1_run1 sample1_run2 sample1_run3 sample2_run1 sample2_run2
      pep_1      17.55388     16.50317     18.70053     15.43822     19.70244
      pep_2      15.83758     24.09699     19.65011     23.27206     18.84700
      pep_3      23.73870     21.87331     17.50231     15.43366     16.56922
      pep_4      18.29231     15.10140     21.22343     24.89333     20.59139
      pep_5      17.22276     23.05993     22.63031     15.89846     16.46224
      pep_6      19.01648     20.47883     19.65011     23.17690     23.13152
      pep_7      15.72499     22.14794     15.71257     24.38314     21.84061
      pep_8      15.02450     21.34943     19.18439     21.20169     18.87477
      pep_9      16.37094     21.18414     18.79415     20.77638     19.04956
      pep_10     16.91909     19.92332     24.08768     20.86509     15.45872
             sample2_run3 sample3_run1 sample3_run2 sample3_run3
      pep_1      17.18089     23.24674     18.46458     18.87797
      pep_2      24.87634     16.16851     20.08279     18.90382
      pep_3      18.48462     21.06622     22.56324     19.38635
      pep_4      19.52186     18.56561     19.84571     19.63025
      pep_5      15.20986     24.65195     18.95437     20.83005
      pep_6      22.49727     24.67561     22.02452     18.87797
      pep_7      16.62958     21.06622     18.33825     15.05981
      pep_8      18.19041     18.32156     17.18413     17.18942
      pep_9      20.97890     21.53844     20.87668     21.00017
      pep_10     20.87997     21.14550     19.84571     18.65062

# BPCA

    Code
      D_BPCA
    Output
                   sample1.pep_1 sample1.pep_2 sample1.pep_3 sample1.pep_4
      sample1_run1      17.41579      15.83758      23.73870      18.29231
      sample1_run2      16.50317      24.09699      21.87331      15.10140
      sample1_run3      18.70053      19.56291      17.50231      21.22343
                   sample1.pep_5 sample1.pep_6 sample1.pep_7 sample1.pep_8
      sample1_run1      17.22276      19.01648      15.72499      15.02450
      sample1_run2      23.05993      20.37257      22.14794      21.34943
      sample1_run3      22.63031      19.56291      15.71257      19.18439
                   sample1.pep_9 sample1.pep_10 sample2.pep_1 sample2.pep_2
      sample1_run1      16.37094       16.91909      15.43822      23.27206
      sample1_run2      21.18414       19.92332      19.70244      18.84700
      sample1_run3      18.79415       24.08768      17.18089      24.87634
                   sample2.pep_3 sample2.pep_4 sample2.pep_5 sample2.pep_6
      sample1_run1      15.43366      24.89333      15.89846      23.17690
      sample1_run2      16.56922      20.59139      16.46224      23.13152
      sample1_run3      18.48462      20.23177      15.20986      22.49727
                   sample2.pep_7 sample2.pep_8 sample2.pep_9 sample2.pep_10
      sample1_run1      24.38314      21.20169      20.68338       20.86509
      sample1_run2      21.84061      18.87477      19.12117       15.45872
      sample1_run3      16.62958      18.19041      20.97890       20.87997
                   sample3.pep_1 sample3.pep_2 sample3.pep_3 sample3.pep_4
      sample1_run1      23.24674      16.16851      20.82845      18.56561
      sample1_run2      18.46458      20.08279      22.56324      19.73177
      sample1_run3      18.73731      18.90382      19.38635      19.63025
                   sample3.pep_5 sample3.pep_6 sample3.pep_7 sample3.pep_8
      sample1_run1      24.65195      24.67561      20.82845      18.32156
      sample1_run2      18.95437      22.02452      18.33825      17.18413
      sample1_run3      20.83005      18.73731      15.05981      17.18942
                   sample3.pep_9 sample3.pep_10
      sample1_run1      21.53844       21.14550
      sample1_run2      20.87668       19.73177
      sample1_run3      21.00017       18.65062

