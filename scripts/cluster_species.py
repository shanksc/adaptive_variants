import argparse
import numpy as np
from sklearn.cluster import KMeans


def get_args():
    parser = argparse.ArgumentParser(description='cluster populations using k means on PC1', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--PCA', type=str, required=True, help='file containing PCA results')
    #this can be expanded to accept multiple samples files which can then support different populations
    parser.add_argument('--n_components', type=int, required=False, help='nth principal component to cluster on', default=1)
    parser.add_argument('--k', type=int, required=True, help='number of clusters')
    parser.add_argument('--o', type=str, required=True, help='output filename')
    
    return parser.parse_args()


def build_X(pcs, n_components):
    
    X = np.zeros((len(pcs), n_components))
    for i in range(len(pcs)):
        for j in range(n_components):
            X[i,j] = pcs[i][j]
    
    return X 


#IID    PC1 
def parse_pca(PCA, n_components):
    
    pcs = []
    labels = []
    with open(PCA, 'r') as f:
        next(f)
        for line in f.readlines():
            split = line.split('\t')
            labels.append(split[0])
            pcs.append([float(split[i+1]) for i in range(n_components)])
    
    return pcs, labels


def cluster(pcs, labels, k, n_components):
    
    X = build_X(pcs, n_components)
    kmeans = KMeans(n_clusters=k, init='k-means++').fit(X)

    return kmeans.labels_


def write_populations(labels, cluster_labels, k, out):
    
    for i in range(k):
        with open(out+f'_{i}.txt', 'w') as f:
            for j in range(len(labels)):
                if cluster_labels[j] == i:
                    f.write(f'{labels[j]}\n')
    
#Could be some benefit to expanding this to more PCs since then we can have more populations
def main():
    
    args = get_args()
    
    #only pcs for now
    pcs, labels = parse_pca(args.PCA, args.n_components)
    #cluster with kmeans
    #print(len(pcs))
    #print(len(labels))
    cluster_labels = cluster(pcs, labels, args.k, args.n_components)
    #print(len(cluster_labels))
    #write a file of the sample names for each population 
    write_populations(labels, cluster_labels, args.k, args.o)
    
    
if __name__ == '__main__':
    main()
