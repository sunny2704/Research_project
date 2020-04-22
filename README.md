# Study of Internal indices on Time Series Data
This repository deals with my research project on Time series.

Clustering is a popular method of unsupervised learning. To find out the goodness of the clustering, we use two types of indices - External Indices and internal indices.
Where external indices have a class label or ground truth to help classify, internal indices rely only on features and quantity inherent to the dataset and they do not have any ground truths with them.
These indices measure the quality or goodness of the clustering based on Euclidean distance but it does not work on Time series data as we use DTW(Dynamic Time Warping distance) for time series.
The purpose of this research is two-fold:
a) To study if we can effectively use DTW instead of Euclidean distance for internal indices.
b) To explore if we have a clear winner between the clustering methods that we are studying.


In the current study, we use 5 cluster algorithms:
1) Kmeans
2) Hierarchical Single Linkage
3) Hierarchical Complete Linkage
4) Hierarchical Average Linkage
5) Distance Density Clustering(DDC) - newly proposed algorithm mentioned in (Ref: Ruizhe Ma and Rafal Angryk, 2017. 2375-9259/17 © 2017 IEEE DOI 10.1109/ICDMW.2017.11)


We study 7 internal indices here - Ball Hall Index, Banfeld Raftery Index, Calinsky Harabasz Index, Davies Bouldin Index, Log SS ratio Index, PBM Index and RTI Index.

The datasets that we have studied are taken from UCR data repository (https://www.cs.ucr.edu/~eamonn/time_series_data_2018/) and they are - ArrowHead, BeetleFly, BirdChicken, Coffee, Gunpoint, Ham, Meat, SonyAIBORobotSurface, ToeSegmentation2 and Wine 
