from assembly import SequenceAssembler

""" Test case #1: small test case from email

    Test case #2: large test case from email

    Test case #3: the first read contains the second read

    ATTAGACCTG
       AGACCT
         ACCTGCC
            TGCCGG
    --------------
    ATTAGACCTGCCGG

    Test case #4: the second read contains the first read

    ATTAGACCTG
    ATTAGACCTGCC
        GACCTGCCTGAC
            TGCCTGACAGC
    -------------------
    ATTAGACCTGCCTGACAGC

    Test case #5: the first and the second reads are equal

    ATTAGACCTG
    ATTAGACCTG
        GACCTGCCTGAC
            TGCCTGACAGC
    -------------------
    ATTAGACCTGCCTGACAGC

    Test case #6: all reads are equal

    ATTAGACCTG
    ATTAGACCTG
    ATTAGACCTG
    ATTAGACCTG
    ----------
    ATTAGACCTG

    Test case #7: the first read can be glued to the right side
                  of the second read

    ATTAGACCAT
        GACCATTAGA
            ATTAGACCGG
                GACCGGATCG
    ----------------------
    ATTAGACCATTAGACCGGATCG

"""

TESTS = [
    ('data/case1_input.txt', 'data/case1_output.txt'),
    ('data/case2_input.txt', 'data/case2_output.txt'),
    ('data/case3_input.txt', 'data/case3_output.txt'),
    ('data/case4_input.txt', 'data/case4_output.txt'),
    ('data/case5_input.txt', 'data/case5_output.txt'),
    ('data/case6_input.txt', 'data/case6_output.txt'),
    ('data/case7_input.txt', 'data/case7_output.txt')
]

def test_case((input_handle, output_handle)):

    s = SequenceAssembler()

    with open(input_handle) as data:
        s.read_fasta(data)

    with open(output_handle) as data:
        output = data.readline().strip()

    s.assemble()

    if s.sequence == output:
        return True
    else:
        return False
    

if __name__ == "__main__":

    passed = 0
    total = len(TESTS)

    for i in range(total):
        if test_case(TESTS[i]):
            passed += 1

    print "Passed {} tests out of {}.".format(passed, total)

