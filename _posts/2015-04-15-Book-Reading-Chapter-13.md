---
layout: post
title: Book Reading -- Chapter 13 Observing and Measuring Social Interaction
date: 2015-04-15
author: Monaen
thumbnail: /images/post_icon/BookReading.jpg
---

## Basic Conception
### Peer effects : 
[If a person's behavior is influenced by one or more other persons and their interactions can be identified with the group effect is there; but here, "others" must be "the same group of persons" (peers), that is to say, these persons are in the same or similar position, and everyone is related with equal relationships.](http://wiki.mbalib.com/wiki/%E5%90%8C%E7%BE%A4%E6%95%88%E5%BA%94)

### Homophily : 
[Isss the tendency of individuals to associate and bond with similar others. The presence of homophily has been discovered in a vast array of network studies.](http://en.wikipedia.org/wiki/Homophily)(i.e., "love of the same") 

### Instrument Variable : 
Instrumental variable methods allow consistent estimation when the explanatory variables (covariates) are correlated with the error terms of a regression relationship. Such correlation may occur when the dependent variable causes at least one of the covariates ("reverse" causation), when there are relevant explanatory variables which are omitted from the model, or when the covariates are subject to measurement error. 

In linear models, there are **two main requirements** for using an IV:

* 1. The instrument must be correlated with the endogenous explanatory variables, conditional on the other covariates.
* 2. The instrument cannot be correlated with the error term in the explanatory equation (conditional on the other covariates), that is, the instrument cannot suffer from the same problem as the original predicting variable.

### Endogeneity : 
In a statistical model, a parameter or variable is said to be endogenous when there is a correlation between the parameter or variable and the error term.

### Community structure : 
In the study of complex networks, a network is said to have community structure if the nodes of the network can be easily grouped into (potentially overlapping) sets of nodes such that each set of nodes is densely connected internally.

iteratively removes edges by calculating the betweenness of each link and then removing the link that has maximum betweenness. The logic is that if a link has a very high betweenness score, then it is connecting (at least) two groups of nodes that otherwise are quite separate.

## Algorithm

### Repeated Bisection : 
one can start with the problem of splitting the set of nodes into some number of (nearly) equal sized groups in a way that minimizes the links across groups and/or maximizes the numb er of links within groups.

### CONCOR Algorithm : 
* 1. leads to a correlation matrix C, where cij is the correlation between row gi and row gj of the adjacency matrix.
* 2. the algorithm then measures the correlation between rows ci and cj of the correlation matrix. In this way, we form a new correlation matrix C(2)
* 3. Iterating
* 4. converges as t grows, so that the entries of C(t) approach some limit.

CONCOR code : 


<pre><code># Stack matricies and matrix transposes
KNOKI <- read.table("KNOKI.csv",header=FALSE,sep=",")
KNOKM <- read.table("KNOKM.csv",header=FALSE,sep=",")
KNOKIT <- t(KNOKI)
KNOKMT <- t(KNOKM)
KNOK <- rbind(KNOKI,KNOKIT,KNOKM,KNOKMT)

# CONCOR function
# CONvergence of iterated CORrelations
# Creates a square matrix of correlations of the column pairs of a matrix
CONCOR <- function(mat){
    colN <- ncol(mat)
    X <- matrix(rep(0,times=colN*colN),nrow=colN,ncol=colN)
    for(i in 1:colN){
        for(j in i:colN){
            X[i,j] <- cor(mat[,i],mat[,j],method=c("pearson"))
        }
    }
    X <- X+(t(X)-diag(diag((X))))
return(X)
}

KNOKSIM <- CONCOR(KNOK)
</code></pre>

The dataset KNOKI can be downloaded here: 

[KNOKI.csv](http://www.casos.cs.cmu.edu/computational_tools/datasets/sets/knokbur/knokbur_organization_organization/knokbur_organization_organization[KNOKI].csv)</br>
[KNOKM.csv](http://www.casos.cs.cmu.edu/computational_tools/datasets/sets/knokbur/knokbur_organization_organization/knokbur_organization_organization[KNOKM].csv)

Results : 
<pre><code>
           [,1]        [,2]        [,3]       [,4]        [,5]        [,6]        [,7]        [,8]        [,9]      [,10]
 [1,] 1.0000000  0.14169568  0.14978617 0.45054945  0.27814338  0.10482848  0.29773258  0.25677630  0.34065934 0.10669059
 [2,] 0.1416957  1.00000000 -0.06131393 0.14169568  0.40350877  0.35043832  0.29730058  0.14306585  0.14169568 0.20680477
 [3,] 0.1497862 -0.06131393  1.00000000 0.04279605 -0.04087596 -0.10206207  0.31622777  0.37500000  0.47075654 0.17108830
 [4,] 0.4505495  0.14169568  0.04279605 1.00000000  0.38310314  0.10482848  0.29773258  0.14978617  0.34065934 0.35772728
 [5,] 0.2781434  0.40350877 -0.04087596 0.38310314  1.00000000  0.31706324  0.32315281 -0.04087596  0.06822385 0.15285570
 [6,] 0.1048285  0.35043832 -0.10206207 0.10482848  0.31706324  1.00000000 -0.08606630  0.06804138 -0.06988566 0.41907904
 [7,] 0.2977326  0.29730058  0.31622777 0.29773258  0.32315281 -0.08606630  1.00000000  0.00000000  0.40599897 0.07728982
 [8,] 0.2567763  0.14306585  0.37500000 0.14978617 -0.04087596  0.06804138  0.00000000  1.00000000  0.25677630 0.29329423
 [9,] 0.3406593  0.14169568  0.47075654 0.34065934  0.06822385 -0.06988566  0.40599897  0.25677630  1.00000000 0.35772728
[10,] 0.1066906  0.20680477  0.17108830 0.35772728  0.15285570  0.41907904  0.07728982  0.29329423  0.35772728 1.00000000
</code></pre>
### [Edge Removal Algorithm:](http://en.wikipedia.org/wiki/Girvan%E2%80%93Newman_algorithm)
iteratively removes edges by calculating the betweenness of each link and then removing the link that has maximum betweenness. The logic is that if a link has a very high betweenness score, then it is connecting (at least) two groups of nodes that otherwise are quite separate.

<pre><code>
#!/usr/bin/env python
import networkx as nx
import math
import csv
import random as rand
import sys

#this method just reads the graph structure from the file
def buildG(G, file_, delimiter_):
    #construct the weighted version of the contact graph from cgraph.dat file
    #reader = csv.reader(open("/home/kazem/Data/UCI/karate.txt"), delimiter=" ")
    reader = csv.reader(open(file_), delimiter=delimiter_)
    for line in reader:
        if len(line) > 2:
            if float(line[2]) != 0.0:
                #line format: u,v,w
                G.add_edge(int(line[0]),int(line[1]),weight=float(line[2]))
        else:
            #line format: u,v
            G.add_edge(int(line[0]),int(line[1]),weight=1.0)

#keep removing edges from Graph until one of the connected components of Graph splits into two
#compute the edge betweenness
def CmtyGirvanNewmanStep(G):
    #print "call CmtyGirvanNewmanStep"
    init_ncomp = nx.number_connected_components(G)    #no of components
    ncomp = init_ncomp
    while ncomp <= init_ncomp:
        bw = nx.edge_betweenness_centrality(G, weight='weight')    #edge betweenness for G
        #find the edge with max centrality
        max_ = max(bw.values())
        #find the edge with the highest centrality and remove all of them if there is more than one!
        for k, v in bw.iteritems():
            if float(v) == max_:
                G.remove_edge(k[0],k[1])    #remove the central edge
        ncomp = nx.number_connected_components(G)    #recalculate the no of components

#compute the modularity of current split
def _GirvanNewmanGetModularity(G, deg_, m_):
    New_A = nx.adj_matrix(G)
    New_deg = {}
    New_deg = UpdateDeg(New_A, G.nodes())
    #Let's compute the Q
    comps = nx.connected_components(G)    #list of components    
    print 'no of comp: %d' % nx.number_connected_components(G)
    Mod = 0    #Modularity of a given partitionning
    for c in comps:
        EWC = 0    #no of edges within a community
        RE = 0    #no of random edges
        for u in c:
            EWC += New_deg[u]
            RE += deg_[u]        #count the probability of a random edge
        Mod += ( float(EWC) - float(RE*RE)/float(2*m_) )
    Mod = Mod/float(2*m_)
    #print "Modularity: %f" % Mod
    return Mod

def UpdateDeg(A, nodes):
    deg_dict = {}
    n = len(nodes)  #len(A) ---> some ppl get issues when trying len() on sparse matrixes!
    B = A.sum(axis = 1)
    for i in range(n):
        deg_dict[nodes[i]] = B[i, 0]
    return deg_dict

#run GirvanNewman algorithm and find the best community split by maximizing modularity measure
def runGirvanNewman(G, Orig_deg, m_):
    #let's find the best split of the graph
    BestQ = 0.0
    Q = 0.0
    while True:    
        CmtyGirvanNewmanStep(G)
        Q = _GirvanNewmanGetModularity(G, Orig_deg, m_);
        print "current modularity: %f" % Q
        if Q > BestQ:
            BestQ = Q
            Bestcomps = nx.connected_components(G)    #Best Split
            print "comps:"
            print Bestcomps
        if G.number_of_edges() == 0:
            break
    if BestQ > 0.0:
        print "Best Q: %f" % BestQ
        print Bestcomps
    else:
        print "Best Q: %f" % BestQ

def main(argv):
    if len(argv) < 2:
        return 1
    graph_fn = argv[1]
    G = nx.Graph()  #let's create the graph first
    buildG(G, graph_fn, ',')

    print G.nodes()
    print G.number_of_nodes()
    
    n = G.number_of_nodes()    #|V|
    A = nx.adj_matrix(G)    #adjacenct matrix

    m_ = 0.0    #the weighted version for number of edges
    for i in range(0,n):
        for j in range(0,n):
            m_ += A[i,j]
    m_ = m_/2.0
    print "m: %f" % m_

    #calculate the weighted degree for each node
    Orig_deg = {}
    Orig_deg = UpdateDeg(A, G.nodes())

    #run Newman alg
    runGirvanNewman(G, Orig_deg, m_)

if __name__ == "__main__":
    sys.exit(main(sys.argv))
</code></pre>

### Hierarchical Clustering : 
One could use the correlation coefficients between the rows of the adjacency matrix, C(1), which was the first step of the CONCOR process above. Instead, we could calculate the
distance between the vectors gi and gj (using Euclidean distance; or city-block distance the number of entries that differ across the two vectors). 

A direct graph illustrate the progress of Hierarchical Clustering Algorithm is as follows:

![Hierarchical Clustering](/images/post_img/chapter13.png)

This picture is measured by similarity between nodes i and j, via a slight variation on the city block distance between rows, or in particular set the distance between nodes i and j to be

![Hierarchical Clustering](/images/post_img/chapter13.1.png)










