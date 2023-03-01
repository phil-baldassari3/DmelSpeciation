import matplotlib.pyplot as plt
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


plt.figure(figsize=(10,8.5))
for idx, i in enumerate(df_ls):

    x = mat_df['x_{id}'.format(id=str(int(idx+1)))].to_list()
    y = mat_df['y_{id}'.format(id=str(int(idx+1)))].to_list()

    if 'm' in i['sex'].to_list():
        s = 'm'
        c = 'blue'
    elif 'f' in i['sex'].to_list():
        s = 'f'
        c = 'red'

    plt.plot(x, y, label = s, color=c)

#plt.legend(bbox_to_anchor=(1, 0.5), loc="center left")

handles, labels = plt.gca().get_legend_handles_labels()
by_label = dict(zip(labels, handles))
plt.legend(by_label.values(), by_label.keys(), bbox_to_anchor=(1, 0.5), loc="center left")

plt.show()


