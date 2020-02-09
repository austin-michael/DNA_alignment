---DNA alignment to study the human family---

Run DNAAlignment.py and a menu will appear letting you know how to use the program.

You have the following options:

- Run unit tests
- Compare Homo Sapien DNA to Neanderthal DNA
- Compare Homo Sapien DNA to Gorilla DNA
- Compare Neanderthal DNA to Gorilla DNA
- Run an algorithm performance test, just input a size and two random DNA sequences will be generated based on that size



FINAL SUMMARY DISCUSSION:

While implementing and testing this algorithm, I have noticed it is extremely efficient. See algorithm_performance.txt for a list
of runtimes (in seconds) given random sequences up to length 900. However, the recursive traceback requires a lot of memory. 
I found sequences with lengths greater than 900 will hit the max recursion depth. I did not, however, attempt to increase my max recursion depth
to test with lengths greater than 900. Based upon this information, I can conclude that the algorithm will only be limited by memory.  