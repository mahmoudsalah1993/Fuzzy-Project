- The theta scopes need to be the sum of each element theta squared, not the vector itself. [[work in diagram]] [[uppdate : done]]
- The simulation works but most graphs are wrong.
- In equation to get mhdot the condition to determine mhdot contains m(error rate) not mh only, error rate = m - mh "below eq 41", and we don't have m as far as i can see
 so error in lines 148 - 153 in code, also 162 - 165.
- The input vector to test with can be found from the sim file, in the integrator block, the initial conditions.
- test from workspace with initial vector : 
ATFC([0.1,0.1,0,0,3.5,3.5,0,0,-0.1,-0.3,-0.5,-0.7,-0.9,0.1,0.3,0.5,0.7,0.9,-0.1,-0.3,-0.5,-0.7,-0.9,0.1,0.3,0.5,0.7,0.9,0])