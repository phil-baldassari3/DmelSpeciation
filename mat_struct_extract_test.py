import scipy.io as spio
import numpy as np
import pandas as pd





mat_dict = spio.loadmat('trx.mat', simplify_cells=True)


#df for each row of struct
df_ls = []
for idx in range(len(mat_dict['trx'])):

    seriesls = []
    for k, v in mat_dict['trx'][idx].items():
        seriesls.append(pd.Series(v, name=k))

    df_ls.append(pd.concat(seriesls, axis=1))


#extract param
param = ['x', 'y']

new_d = {}

for p in param:
    for i in df_ls:
        l = i[p].to_list()
        new_d.update({p + '_' + str(int((i['id'].to_list()[0]))) : l})



mat_df = pd.DataFrame(new_d)
mat_df.to_csv("xy_mat_test.csv", index=False)











perframe_dict = spio.loadmat('perframe/dist_to_wall.mat', simplify_cells=True)


new_perframe = {}
for idx in range(len(perframe_dict['data'])):
    new_perframe.update({idx+1 : perframe_dict['data'][idx]})


perframe_df = pd.DataFrame(new_perframe)
perframe_df.to_csv("perframe_dist2wall_test.csv", index=False)