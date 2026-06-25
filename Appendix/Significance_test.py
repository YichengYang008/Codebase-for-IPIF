import numpy as np
import scipy.stats as stats
from statsmodels.stats.multitest import multipletests
import pandas as pd
import seaborn as sns

# 示例 RMSE 数据
rmse_matrix = np.array([
   [0.0539451, 0.0543201, 0.0539899, 0.0539971, 0.0536891],  
   [0.0630414, 0.0640897, 0.063841, 0.0637934, 0.0634199],  
   [0.27302, 0.270247, 0.268085, 0.275564, 0.270798],  
  [0.063142, 0.0641825, 0.0639258, 0.0638741, 0.0635129]   
])

method_names = ['UP-FHDI', 'NAIVE', 'GAIN', 'HI-VAE']
n_methods = len(method_names)

print(rmse_matrix[0])
# 1. Friedman 检验
friedman_stat, friedman_p = stats.friedmanchisquare(*rmse_matrix)
print(f"\nFriedman Test p-value: {friedman_p:.4f}")

# 2. 若显著则执行 A vs B/C/D 的 Wilcoxon 单尾检验
if friedman_p < 0.05:
    print("\nFriedman 检验显著，进行 UP-FHDI vs 其他方法的 Wilcoxon 单尾检验（Bonferroni 校正）...\n")
    comparisons = []
    p_vals = []


    for j in range(1, n_methods):  # A vs B/C/D
        
        w_stat, p_val = stats.wilcoxon(rmse_matrix[0], rmse_matrix[j], alternative='less')
        #w_stat, p_val = stats.ttest_rel(rmse_matrix[0], rmse_matrix[j], alternative='less')
        print(p_val)
        comparisons.append({
            'Method_A': 'UP-FHDI',
            'Method_Other': method_names[j],
            'Raw p-value': p_val
        })
        p_vals.append(p_val)
    
    df_comp1 = pd.DataFrame(comparisons)
    print(df_comp1)
    #df_comp1.to_csv("ttest1.csv", index=False)


    # 多重比较校正（Bonferroni）
    reject, corrected_pvals, _, _ = multipletests(p_vals, method='bonferroni')
    for idx, comp in enumerate(comparisons):
        comp['Bonferroni-corrected'] = corrected_pvals[idx]
        comp['Significant (Bonferroni < 0.05)'] = reject[idx]

    df_comp2 = pd.DataFrame(comparisons)
    print(df_comp2)


    #df_comp2.to_csv("ttest2.csv", index=False)

else:
    print("\nFriedman 检验未显著，不做 Wilcoxon 检验。")



#情况	                      正态？	使用方法
#同一数据、不同方法（paired）	✅ 是	Paired t-test
#同一数据、不同方法（paired）	❌ 否	Wilcoxon signed-rank test
