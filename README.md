Graph Coarsening

- Serial version of Graph Coarsening Algorithm proposed by: "Parallel, Portable Algorithms for Distance-2
Maximal Independent Set and Graph Coarsening" -- Brian Kelley and Sivasankaran Rajamanickam -- Sandia National Laboratories, Albuquerque, New Mexico, U.S.A

- csv = csv files of graphs
- suitesparse = suitesparse graphs downloaded
- convert_graph_formats.py = reads in files to graphs and converts different data structures
- coarse_serial.py = serial version of coarsening algorithm proposed in paper
    - given an input graph and a "goal graph size", it will output a coarsened graph that is of size <= to the given goal size


- graph_coase_simple.py = simple, serial graph coarsening algorithm based on "maximal matching" dicussed in "A Multilevel Algoritm for Partitioning Graphs" -- BRUCE HENDRICKSON AND ROBERT LELAND, SANDIA NATIONAL LABORATORIES

** Needs verification mechanisms **
