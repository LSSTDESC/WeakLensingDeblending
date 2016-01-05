#return sum of n, n-1, ..., k
def sum_n_to_k(n,k):
    l = [x for x in range(n+1) if x>=k]
    return sum(l)

"""
these function returns the correction position in the datacube of the specified partial derivative
i,j. Note that we are filling in cubes number 7 to 20 and partials commute. Order established in 
the variations dict.
Differentiating partial flux with respect to flux gives 0. 
Differentiating any partial with respect to flux leaves it unchanged. 
So we do not have to calculate any partials with respect to flux at all. 
"""
def second_partials_indices_to_datacubes_indices(i,j):
    return sum_n_to_k(5,7-i) + j - i + 7 


# Prepare the datacube that we will return.

if no_partials and no_bias:
    ncube = 1

elif not no_partials and no_bias:
    # The nominal image doubles as the flux partial derivative.
    ncube = 1+len(variations)

else:
    ncube = 1 + len(variations) + 15 #number of second partials for 5 parameters (no flux).

height,width = cropped_stamp.array.shape
datacube = np.empty((ncube,height,width))
datacube[0] = cropped_stamp.array #flux partial is the same. 

#calculate partials and second partials, if requested. 
if not no_partials and not no_bias:
    for i,(pname_i,delta_i) in enumerate(variations):

        #partial
        variation_stamp = (galaxy.renderer.draw(**{pname: +delta_i}).copy() - 
                               galaxy.renderer.draw(**{pname: -delta_j}))
        datacube[i+1] = variation_stamp.array/(2*delta_i)

        if not no_bias:
        
            for j,(pname_j,delta_j) in range(len(variations))[i:],variations:

                ##why the .copy()?? 
                galaxy_iup_jup = galaxy.renderer.draw(**{pname_i: +delta_i, pname_j: +delta_j})
                galaxy_iup_jdown = galaxy.renderer.draw(**{pname_i: +delta_i, pname_j: -delta_j})
                galaxy_idown_jup = galaxy.renderer.draw(**{pname_i: -delta_i, pname_j: +delta_j})
                galaxy_idown_jdown = galaxy.renderer.draw(**{pname_i: -delta_i, pname_j: -delta_j})

                partial_i_j = galaxy_iup_jup + galaxy_idown_jdown - galaxy_idown_jup - galaxy_iup_jdown
                datacube[f(i,j)] = partial_i_j

# Calculate partial derivative images only, if requested.
if not no_partials:
    for i,(pname,delta) in enumerate(variations):
        variation_stamp = (galaxy.renderer.draw(**{pname: +delta}).copy() - 
            galaxy.renderer.draw(**{pname: -delta}))
        datacube[i+1] = variation_stamp.array/(2*delta)
