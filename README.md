### 1. Building the Artifact


We need a linux operating system with build-essential and cmake installed. Run the following commands to build:
```shell script
$ git clone git@github.com:qingkaishi/context-sensitive-reachability.git csr
$ mkdir csr/build && cd csr/build
$ cmake .. -DCMAKE_BUILD_TYPE=Release
$ make
```

The commands above produces an executable `csr` to test how existing indexing schemes (pathtree and grail) speed up
context-sensitive reachability. Read the following paper for more information on indexing context-sensitive reachability:
```
OOPSLA'21: Indexing Context-Sensitive Reachability
Qingkai Shi, Yongchao Wang, Charles Zhang
The 36th ACM SIGPLAN Conference on Objected Oriented Programming, Systems, Languages, and Applications
```

### 2. Dataset

We include some data in the repo (see the folder `dataset/`) for testing our artifact.
Each file in the folder `dataset/` represents a program-valid graph. Each graph file is in the following format, e.g.,
```text
graph_for_greach
6
0: 1 2.1 #0
1: 2.2 #0
2: 3 #1
3: 4.-1 5.-2 #1
4: #0
5: #0
```

The first line of the graph file is just a string `graph_for_greach` to indicate the format of the file.

The second line indicates the number of vertices in the graph.

The following lines represents the graph as an adjacency list where the vertices are consecutive natural numbers,
which are in the format below: 
```text
source: (target[.[-]call_id] )*#function_id
```
For example:
`0: 1 2.1 #0` means the vertex `0` is in the function `0` and there are two edges:
* the edge from `0` to `1`, which is an intra-procedural edge
* the edge from `0` to `2`, which is a call edge at the call site `1`.

For example:
`3: 4.-1 5.-2 #1` means the vertex `3` is in the function `1` and there are two edges:
* the edge from `3` to `4`, which is a return edge at the call site `1`.
* the edge from `3` to `5`, which is a return edge at the call site `2`.

### 3. Running the Experiments


```
$ ./csr -h

Usage:
        csr [-h] [-t] [-m pathtree_or_grail] [-d grail_dim] [-n num_query] [-q query_file] [-g query_file] graph_file
Description:
        -h      Print the help message.
        -n      # reachable queries and # unreachable queries to be generated, 100 for each by default.
        -g      Save the randomly generated queries into file.
        -q      Read the randomly generated queries from file.
        -t      Evaluate transitive closure.
        -r      Evaluate rep's tabulation algorithm.
        -m      Evaluate what indexing approach, pathtree, grail, or pathtree+grail.
        -d      Set the dim of Grail, 2 by default.
```

Sample usage:

```
$ ./csr -n 1000 -m grail ./dataset/mcf.txt

.....
.....
.....
--------- GRAIL Queries Test ------------
Grail for 1000 reachable queries: 4 ms. Success rate: 100 %.
Grail for 1000 unreachable queries: 0 ms. Success rate: 100 %.
--------- Indexing Construction Summary ---------
#Vertices: 22194
#Edges: 29372
#Summary Edges: 2452
Summary Edge     time: 10.58 ms. 
Summary Edge     size: 0.02 mb.
GRAIL    indices time: 36.76 ms. 
GRAIL    indices size: 1.19 mb. 
```

### 4. Acknowledgement

This repo includes the source code contributed by the authors of [PathTree](http://www.cs.kent.edu/~nruan/soft.html) and [Grail](https://github.com/zakimjz/grail). 
We deeply appreciate them contributing to the code.