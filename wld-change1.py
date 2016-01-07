
"""
This function returns the correction position in the datacube of the specified 2nd 
partial derivative i,j. Note that we are filling in cubes index 6 to 19 
Follow the order in the variations dict. 
Differentiating partial flux with respect to flux gives 0. 
Differentiating any partial with respect to flux leaves it unchanged. 
Therefore ignore flux partials.
"""
def index_partial(i,j, variations):
    
    #return sum of n, n-1, ..., k (non-excluding)
    def sum_n_to_k(n,k):
        l = [x for x in range(n+1) if x>=k]
        return sum(l)

    return sum_n_to_k(len(variations),len(variations)+2-i) + j - i + (1+len(variations)) 

# Prepare the datacube that we will return.
ncube = 1

if no_partials:
    ncube = 1+len(variations)
    if no_bias:
         #15 is number of second partials for 5 parameters (no flux).
        ncube = 1 + len(variations) + 15


height,width = cropped_stamp.array.shape
datacube = np.empty((ncube,height,width))
datacube[0] = cropped_stamp.array #flux partial is the same. 

#calculate partials, if requested.
if not no_partials:
      #The nominal image doubles as the flux partial derivative.
    for i,(pname_i,delta_i) in enumerate(variations):

        variation_stamp = (galaxy.renderer.draw(**{pname: +delta_i}).copy() - 
                               galaxy.renderer.draw(**{pname: -delta_i}))
        datacube[i+1] = variation_stamp.array/(2*delta_i)

        #calculated second partials, if requested. 
        if not no_bias:        
            for j,(pname_j,delta_j) in range(len(variations))[i:],variations:
                galaxy_iup_jup = galaxy.renderer.draw(**{pname_i: +delta_i, pname_j: +delta_j})
                galaxy_iup_jdown = galaxy.renderer.draw(**{pname_i: +delta_i, pname_j: -delta_j})
                galaxy_idown_jup = galaxy.renderer.draw(**{pname_i: -delta_i, pname_j: +delta_j})
                galaxy_idown_jdown = galaxy.renderer.draw(**{pname_i: -delta_i, pname_j: -delta_j})

                variation_i_j = (galaxy_iup_jup + galaxy_idown_jdown - 
                                 galaxy_idown_jup - galaxy_iup_jdown)
                datacube[index_partial(i+1,j+1, variations)] = variation_i_j / (4*delta_i*delta_j)