import numpy as np
import random
import time as t

# global lists used to track the alignment from the recursive traceback. (.clear() is called before every usage)
alignmentA = []
alignmentB = []


# A = Sequence 1
# B = Sequence 2
# n = length of A
# m = length of B
# edit_dict is dictionary of scores to assign
def dna_alignment(n, m, A, B, edit_dict):
    cache = np.zeros((n + 1, m + 1))

    # fill in the base cases (i.e. B is empty or A is empty)
    for i in range(1, n):
        cache[i][0] = cache[i - 1][0] + assign_edit_distance(A[i], '-', edit_dict)

    for j in range(1, m):
        cache[0][j] = cache[0][j - 1] + assign_edit_distance('-', B[j], edit_dict)

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            # If the last characters are different, edit is needed, test inserting, removing, and replacing. Return the max and store it in the cache
            cache[i][j] = max(cache[i][j - 1] + assign_edit_distance('-', B[j - 1], edit_dict),  # insert
                              cache[i - 1][j] + assign_edit_distance(A[i - 1], '-', edit_dict),  # remove
                              cache[i - 1][j - 1] + assign_edit_distance(A[i - 1], B[j - 1], edit_dict))  # replace

    # These four variable declarations are unnecessary, however, they aid in understanding the traceback structure.
    i = n
    j = m
    A_i = n - 1
    B_j = m - 1
    traceback(i, j, A_i, B_j, cache, A, B, edit_dict)
    return int(cache[n][m])


def traceback(i, j, A_i, B_j, cache, A, B, edit_dict):

    if i == 0:  # sequence A is empty, finish B
        for k in range(0, j):
            alignmentA.append('-')
            alignmentB.append(B[k])
        return
    elif j == 0:  # sequence B is empty, finish A
        for k in range(0, i):
            alignmentA.append(A[k])
            alignmentB.append('-')
        return
    elif cache[i][j] == cache[i][j - 1] + assign_edit_distance(A[A_i], B[B_j - 1], edit_dict):  # test if it was an insertion
        alignmentA.append(A[A_i])
        alignmentB.append(B[B_j - 1])
        return traceback(i, j - 1, A_i, B_j - 1, cache, A, B, edit_dict)
    elif cache[i][j] == cache[i - 1][j] + assign_edit_distance(A[A_i - 1], B[B_j], edit_dict):  # test if it was a removal
        alignmentA.append(A[A_i - 1])
        alignmentB.append(B[B_j])
        return traceback(i - 1, j, A_i - 1, B_j, cache, A, B, edit_dict)
    else: #cache[i][j] == cache[i - 1][j - 1] + assignEditDistance(A[A_i - 1], B[B_j - 1], edit_dict):
        # if it was not an insertion or removal, it was a replacement
        alignmentA.append(A[A_i - 1])
        alignmentB.append(B[B_j - 1])
        return traceback(i - 1, j - 1, A_i - 1, B_j - 1, cache, A, B, edit_dict)


def create_table():

    edit_dict = {
        'AA': 5,
        'AC': -1,
        'AG': -2,
        'AT': -1,
        'A-': -3,
        'CC': 5,
        'CG': -3,
        'CT': -2,
        'C-': -4,
        'GG': 5,
        'GT': -2,
        'G-': -2,
        'TT': 5,
        'T-': -1,
        '--': 0
    }

    return edit_dict


# searches the edit dictionary to assign the corresponding edit distance
# edit_dict does not contain same value pairs (i.e. AC and CA are equivalent)
# first checks the dictionary if the value exists, if it does not it swaps the characters and checks if that pair exists.
# returns 0 if data has bad characters
def assign_edit_distance(A_char, B_char, edit_dict):

    comparison = str(A_char) + str(B_char)
    swapped = str(B_char) + str(A_char)

    if comparison in edit_dict:
        return edit_dict[comparison]
    elif swapped in edit_dict:
        return edit_dict[swapped]
    else:
        return 0


def unit_tests():
    edit_dict = create_table()

    print("_______________________________UNIT TESTS_______________________________")

    A = 'AAAAA'
    B = 'AAAA'
    n = A.__len__()
    m = B.__len__()
    alignmentA.clear()
    alignmentB.clear()
    print(" The expected number of edit steps for the following DNA sequences is 17\n",
          "A: ", A, "\n",
          "B: ", B, "\n",
          "Score: ", dna_alignment(n, m, A, B, edit_dict), "\n")

    print('Alignment A: ', alignmentA)
    print('Alignment B: ', alignmentB)

    print("________________________________________________________________________")

    A = 'AAAAA'
    B = 'AAAAG'
    n = A.__len__()
    m = B.__len__()
    alignmentA.clear()
    alignmentB.clear()
    print(" The expected number of edit steps for the following DNA sequences is 18\n",
          "A: ", A, "\n",
          "B: ", B, "\n",
          "Score: ", dna_alignment(n, m, A, B, edit_dict), "\n")

    print('Alignment A: ', alignmentA)
    print('Alignment B: ', alignmentB)

    print("________________________________________________________________________")

    A = 'ACGTACGT'
    B = 'ACGT'
    n = A.__len__()
    m = B.__len__()
    alignmentA.clear()
    alignmentB.clear()
    print(" The expected number of edit steps for the following DNA sequences is 10\n",
          "A: ", A, "\n",
          "B: ", B, "\n",
          "Score: ", dna_alignment(n, m, A, B, edit_dict), "\n")

    print('Alignment A: ', alignmentA)
    print('Alignment B: ', alignmentB)
    print("________________________________________________________________________")

    A = 'AAAAACG'
    B = 'GAAAAAG'
    n = A.__len__()
    m = B.__len__()
    alignmentA.clear()
    alignmentB.clear()
    print(" The expected number of edit steps for the following DNA sequences is 23\n",
          "A: ", A, "\n",
          "B: ", B, "\n",
          "Score: ", dna_alignment(n, m, A, B, edit_dict), "\n")

    print('Alignment A: ', alignmentA)
    print('Alignment B: ', alignmentB)
    print("________________________________________________________________________")

    A = 'AAAGCT'
    B = 'CGTACG'
    n = A.__len__()
    m = B.__len__()
    alignmentA.clear()
    alignmentB.clear()
    print(" The expected number of edit steps for the following DNA sequences is 2\n",
          "A: ", A, "\n",
          "B: ", B, "\n",
          "Score: ", dna_alignment(n, m, A, B, edit_dict), "\n")

    print('Alignment A: ', alignmentA)
    print('Alignment B: ', alignmentB)
    print("________________________________________________________________________")

    A = 'AAAGCTTTTT'
    B = 'CGTACG'
    n = A.__len__()
    m = B.__len__()
    alignmentA.clear()
    alignmentB.clear()
    print(" The expected number of edit steps for the following DNA sequences is -2\n",
          "A: ", A, "\n",
          "B: ", B, "\n",
          "Score: ", dna_alignment(n, m, A, B, edit_dict), "\n")

    print('Alignment A: ', alignmentA)
    print('Alignment B: ', alignmentB)
    print("________________________________________________________________________")


def clean_data(data):
    for i in data:
        data = i.split()

    data = [x for x in data if not any(c.isdigit() for c in x)]     # Removes numbers
    data = [x.upper() for x in data]    # converts all lowercase letters to uppercase letters
    del data[0]     # Removes introductory string (i.e. "ORIGIN")
    return data


def compare_dna_sequences(file_1, file_2):
    DNA_A = []
    DNA_B = []
    score = []
    edit_dict = create_table()
    alignmentA.clear()
    alignmentB.clear()

    f = open(file_1, "r")
    DNA_A.append(f.read())
    f.close()

    DNA_A = clean_data(DNA_A)

    f = open(file_2, "r")
    DNA_B.append(f.read())
    f.close()

    DNA_B = clean_data(DNA_B)

    # check A > B or B > A, if true, fill in A or B with empty characters
    if DNA_A.__len__() > DNA_B.__len__():
        for i in range(DNA_B.__len__(), DNA_A.__len__()):
            DNA_B.append('----------')
    if DNA_B.__len__() > DNA_A.__len__():
        for i in range(DNA_A.__len__(), DNA_B.__len__()):
            DNA_A.append('----------')

    # run the algorithm over all sequences. A and B have the same size, so we can just run it over length of A
    for i in range(0, DNA_A.__len__()):
        score.append(dna_alignment(DNA_A[i].__len__(), DNA_B[i].__len__(), DNA_A[i], DNA_B[i], edit_dict))

    # generate total score based upon all scores returned for the sequences
    total_score = 0
    for i in score:
        total_score += i

    # write output to file
    file_name = 'comparing_' + file_1.split('.', 1)[0] + '_' + file_2.split('.', 1)[0] + '.txt'
    f = open(file_name, "w")
    f.write('Score: ' + str(total_score) + '\n')
    f.write('Alignment:' + '\n')
    for i in range(0, DNA_A.__len__()):
        f.write('\t\t\t' + str(alignmentA[i]) + ' --> ' + str(alignmentB[i]) + '\n')
    f.close()


def algorithm_performance_test(length):
    nucleobases = ['A', 'C', 'G', 'T']
    DNA_A = []
    DNA_B = []
    alignmentA.clear()
    alignmentB.clear()

    for i in range(length):
        rand_1 = random.randint(0, 3)
        rand_2 = random.randint(0, 3)
        DNA_A.append(nucleobases[rand_1])
        DNA_B.append(nucleobases[rand_2])

    A = ''.join(DNA_A)
    B = ''.join(DNA_B)
    n = A.__len__()
    m = B.__len__()
    edit_dict = create_table()

    begin_time = t.time()
    print("\n A: ", A, "\n",
          "B: ", B, "\n",
          "Score: ", dna_alignment(n, m, A, B, edit_dict), "\n")
    end_time = t.time()
    total_time = end_time - begin_time

    print(' Alignment A: ', alignmentA)
    print(' Alignment B: ', alignmentB)
    print(' Runtime in seconds: ', total_time)


    f = open('algorithm_performance.txt', "a+")
    f.write('Runtime given sequences of length: ' + str(length) + ' was: ' + str(total_time) + '\n')
    f.close()


def run(startup):
    if startup:
        print("_______________________________________________________________________\n")
        print("Welcome to my implementation of DNA alignment to study the human family")
        print("_______________________________________________________________________")

    selection = input("\nEnter 0 for unit tests\n"
                      "Enter 1 to compare Homo sapien DNA to Neanderthal DNA\n"
                      "Enter 2 to compare Homo sapien DNA to Gorilla DNA\n"
                      "Enter 3 to compare Neanderthal DNA to Gorilla DNA\n"
                      "Enter 4 for algorithm performance test\n"
                      "Enter E to exit program\n\n"
                      "NOTE: Output in .txt files will not generate until the program exits.\n\n")

    if selection == '0':
        unit_tests()
        run(False)
    elif selection == '1':
        compare_dna_sequences('Homosapien.txt', 'Neanderthal.txt')
        run(False)
    elif selection == '2':
        compare_dna_sequences('Homosapien.txt', 'Gorilla.txt')
        run(False)
    elif selection == '3':
        compare_dna_sequences('Neanderthal.txt', 'Gorilla.txt')
        run(False)
    elif selection == '4':
        length = input('Enter how long you want the 2 random DNA sequences to be: \n')
        if int(length) >= 1:
            algorithm_performance_test(int(length))
        else:
            print('Length must be greater than 1\n')
        run(False)
    elif selection == 'E' or selection == 'e':
        return
    else:
        print(selection, 'is an invalid character:')
        run(False)


run(True)
