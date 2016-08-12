BEGIN {
#	z_bin_size = 4.0e0
#	z_start = -7.5e0
#The bin size in the y-direction 
#is the entire box length.
	#x_bin_index = 28
#z_bin_index is passed depending on gap width.
#y_bin_index is 1.
	count = 0
	for(i=1;i<=x_bin_index;i++)
	{
		x[i] = x_start + x_bin_size*(i-1)
		density[i] = 0.0
		px[i] = 0.0
		binvol[i] = x_bin_size*y_bin_size*(10.0 + 7.05e-4*(x[i]-300.e0)^2)
		binvol[i] = x_bin_size*y_bin_size*12.0
	}
	#printf("Starting...")	
}
{
	

	if(NF==3)
			count = count + 1

	if(NF==9)
	{
			
		x_level = int(($2 - x_start)/x_bin_size) + 1
		#if(!($4 - h[i] - 24.84 <= 2.e0 || $4 - 22.84 <= 2.e0))
		if($4 >= 27 && $4 <= 35)
		{
				density[x_level] = density[x_level] + $5
				vx[x_level] = vx[x_level] + $5*$6
		}			

			
	}
	

}
END {
	for(i=1;i<=x_bin_index;i++)
	{
			if(density[i] -0.0 > 1.e-5)
			{
					vx[i] = vx[i]/count
					density[i] = density[i]/count
					vx[i] = vx[i]/binvol[i]
			}

			printf("%f %f %f\n", x[i], \
 			density[i], vx[i]) 
	}


}
