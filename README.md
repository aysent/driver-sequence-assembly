# driver-sequence-assembly
Assembly of chromosome sequence from shotgun sequencing fragments

Problem of assembling the original chromosome (sequence) from
its multiple fragments (reads) is represented with a graph,
where vertices are individual reads and edges are overlaps
between reads. The assembly of the original sequence
is equivalent to finding such a path through the graph
that each read is only used once. The path is found using
the Depth First Search (DFS) algorithm.


An example of usage of this module:

    from assembly import SequenceAssembler
    
    s = SequenceAssembler()
    
    # read data
    with open("data.fasta") as data:
        s.read_fasta(data)
    
    # run assembly
    s.assemble()
    
    # access path and resulting sequence
    
    print s.path
    print s.sequence

Run tests:

    python test.py
