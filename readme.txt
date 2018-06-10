- The theta scopes need to be the sum of each element theta squared, not the vector itself. [[work in diagram]]
- The simulation works but most graphs are wrong.
- In equation to get mhdot the condition to determine mhdot contains m(error rate) not mh only, error rate = m - mh "below eq 41", and we don't have m as far as i can see
 so error in lines 148 - 153 in code, also 162 - 165.
- The input vector to test with can be found from the sim file, in the integrator block, the initial conditions.