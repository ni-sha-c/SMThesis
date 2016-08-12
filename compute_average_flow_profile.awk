BEGIN {
#	z_bin_size = 5.e0
#	z_start = -7.5e0
#The bin size in the y-direction 
#is the entire box length.
	#x_bin_index = 28
#z_bin_index is passed depending on gap width.
#y_bin_index is 1.
	count = 0
	for(i=1;i<=z_bin_index;i++)
	{
		z[i] = z_start + z_bin_size*(i-1)
		density[i] = 0.0
		vx[i] = 0.0
	}
	
}
{
	if(NF==2)
			count = count + 1

	if(NF==7)
	{
		
		z_level = int(($1 - z_start)/z_bin_size) + 1
		density[z_level] = density[z_level] + $6
		vx[z_level] = vx[z_level] + $7			

			
	}

}
END {
	for(i=1;i<=z_bin_index;i++)
	{
			if(density[i] -0.0 > 1.0)
			{
					vx[i] = vx[i]/count
					density[i] = density[i]/count
					#vx[i] = vx[i]/density[i]
				printf("%f %f %f\n", z[i], \
 			density[i], vx[i]) 

			}

				}


}
