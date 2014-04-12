SIMPLE
======

Statistical Inference with Multi PLanet Exemplars
Version 0.1.1

simple_lib : Library Module for SIMPLE




Performance tests as of 04/11/14

Use this command
python ABC.py stars.csv simplest_model test_cat.dat 100

test_cat.dat too big for repo, contact me for copy

original code (master branch)

    100 trials of 952 trials accepted
    Total: 1234.82149792s  Catalog Gen + ABC: 1214.74370599s
    Overhead: 20.0777919292s  Time/Catalog: 1.2759912878s
    Total Time: 1236.38190913




Baseline after converted to 1 processor map

    100 trials of 1000 trials accepted
    Total: 1050.59805489s  Catalog Gen + ABC: 1030.92166185s
    Overhead: 19.6763930321s  Time/Catalog: 1.03092166185s
    Total Time: 1051.61795092




parallel map '4 cores' on my machine

    100 trials of 1000 trials accepted
    Total: 656.568832874s  Catalog Gen + ABC: 636.436167955s
    Overhead: 20.1326649189s  Time/Catalog: 0.636436167955s
    Total Time: 657.797276974



On rcc cluster:

    Ran on 16 Processors
    116 trials of 1000 trials accepted
    Total: 287.09539485s  Catalog Gen + ABC: 263.155495882s
    Overhead: 23.9398989677s  Time/Catalog: 0.263155495882s
    Total Time: 287.91439414

