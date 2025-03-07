#importing modules
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score

#kmeans clustering function
def kmeans_data_out(xy, k):
    """
    Inputs:
    xy: (array-like) 2D array of data points from make_blobs()
    k: (int) number of clusters for k-means clustering
    returns list of categories
    """

    #defining clusters
    kmeans = KMeans(init="random", n_clusters=k, n_init=20, max_iter=300)

    #clustering data
    km_prediction = kmeans.fit_predict(xy)

    #sse
    sse_score = kmeans.inertia_

    #silhouette score
    sil_score = silhouette_score(xy, km_prediction)

    return km_prediction, sse_score, sil_score


#SET THESE
input_csv = "PCA_maf0.05_Autosomes.csv"


#opening data
df = pd.read_csv(input_csv)
xy_coords = df[["EV1","EV2"]].to_numpy()

"""
#looping through ks
sse_ls = []
sil_ls = []
for k in range(2, 22):
    clusters, sse, sil = kmeans_data_out(xy_coords, k)
    sse_ls.append(sse)
    sil_ls.append(sil)

plt.plot(range(2, 22), sse_ls)
plt.xticks(range(2, 22))
plt.xlabel("k")
plt.ylabel("WSS")
plt.title("Within-Cluster Sum of Squares")
#plt.show()
plt.savefig(input_csv.replace(".csv", "_elbow.png"))

plt.clf()

plt.plot(range(2, 22), sil_ls)
plt.xticks(range(2, 22))
plt.xlabel("k")
plt.ylabel("Silhouette Score")
plt.title("Silhouette Method")
#plt.show()
plt.savefig(input_csv.replace(".csv", "_sil.png"))





"""
#generating new csv

for final_k in range(3,6):

    clusters, sse, sil_score = kmeans_data_out(xy_coords, final_k)

    df["cluster"] = clusters



    output_csv = input_csv.replace(".csv", "_k{}.csv".format(final_k))
    df.to_csv(output_csv, index=False)

