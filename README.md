
This repository gives the software implementation of the improved 4-way and unbalanced 5-way polynomial multiplication algorithms that are described in the following paper by Yeniaras and Cenk: https://link.springer.com/chapter/10.1007/978-3-031-17510-7_9

The improved algorithms are applied to post-quantum, NIST PQC candidate NTRU Prime KEM by Bernstein et al. and the C implementation provide 29.85%  speedup for sntrup653 and a 35.52% speedup for sntrup761. Furtermore we observe 48.6%  decrease in arithmetic complexity for polynomial multiplication over F_9 and a 26.8% decrease for polynomial multiplication over F_3 respectively.  Please refer to the paper for the details of the algorithms and the hybrid methods we used in the application to NTRU Prime KEM.
