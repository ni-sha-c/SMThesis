/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Nisha Chandramoorthy, MIT (email:nishac at mit dot edu)
------------------------------------------------------------------------- */

#include <mpi.h>
#include <string.h>
#include "compute_stress_mop.h"
#include "atom.h"
#include "update.h"
#include "force.h"
#include "pair.h"
#include "neighbor.h"
#include "neigh_request.h"
#include "neigh_list.h"
#include "group.h"
#include "kspace.h"
#include "error.h"
#include <math.h>
#include "comm.h"
#include "domain.h"
#include "math_const.h"
#include <iostream>
#include <fstream>
using namespace LAMMPS_NS;
using namespace MathConst;

#define SMALL 0.00001

/* ---------------------------------------------------------------------- */

ComputeStressMOP::ComputeStressMOP(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 5) error->all(FLERR,"Illegal compute stress/mop command");

  scalar_flag = 1;
	vector_flag = 1;
  extscalar = 1;
  extvector = 1;
	
  

  zpoint = force->numeric(FLERR,arg[3]);
  xbin_size = force->numeric(FLERR,arg[4]);
 
  pairflag = 1;
  
  int iarg = 4;
  /*
  while (iarg < narg) {
    if (strcmp(arg[iarg],"pair") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal compute group/group command");
      if (strcmp(arg[iarg+1],"yes") == 0) pairflag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) pairflag = 0;
      else error->all(FLERR,"Illegal compute group/group command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"kspace") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal compute group/group command");
      if (strcmp(arg[iarg+1],"yes") == 0) kspaceflag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) kspaceflag = 0;
      else error->all(FLERR,"Illegal compute group/group command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"boundary") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal compute group/group command");
      if (strcmp(arg[iarg+1],"yes") == 0) boundaryflag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) boundaryflag  = 0;
      else error->all(FLERR,"Illegal compute group/group command");
      iarg += 2;
    } else error->all(FLERR,"Illegal compute group/group command");
  }
  */
 
	size_vector = ceil(domain->boxhi[0]/xbin_size); 
	maxsize=size_vector;
  vector = new double[maxsize];
}

/* ---------------------------------------------------------------------- */

ComputeStressMOP::~ComputeStressMOP()
{
  delete [] vector;
}

/* ---------------------------------------------------------------------- */

void ComputeStressMOP::init()
{
  // if non-hybrid, then error if single_enable = 0
  // if hybrid, let hybrid determine if sub-style sets single_enable = 0

  if (pairflag && force->pair == NULL)
    error->all(FLERR,"No pair style defined for compute group/group");
  if (force->pair_match("hybrid",0) == NULL && force->pair->single_enable == 0)
    error->all(FLERR,"Pair style does not support compute group/group");

   if (pairflag) {
    pair = force->pair;
    cutsq = force->pair->cutsq;
  } else pair = NULL;

  if (pairflag) {
    int irequest = neighbor->request(this,instance_me);
    neighbor->requests[irequest]->pair = 0;
    neighbor->requests[irequest]->compute = 1;
    neighbor->requests[irequest]->occasional = 1;
  }
}

/* ---------------------------------------------------------------------- */

void ComputeStressMOP::init_list(int id, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

double ComputeStressMOP::compute_scalar()
{
  invoked_scalar = invoked_vector = update->ntimestep;

  scalar = 0.0;
 	for(int kk=0; kk < maxsize; kk++)
		vector[kk] = 0.e0; 

  pair_contribution();
  

  return scalar;
}

/* ---------------------------------------------------------------------- */

void ComputeStressMOP::compute_vector()
{
  invoked_scalar = invoked_vector = update->ntimestep;

  scalar = 0.0;
  for(int ii=0; ii < maxsize; ii++)
		vector[ii] = 0.e0;

  pair_contribution();
  //if (kspaceflag) kspace_contribution();
}

/* ---------------------------------------------------------------------- */

void ComputeStressMOP::pair_contribution()
{
  int i,j,ii,jj,kk,inum,jnum,itype,jtype,sgn_i,sgn_j;
  double xtmp_i,ytmp_i,ztmp_i,xtmp_j,ytmp_j,ztmp_j,delx,dely,delz;
  double rsq,eng,fpair,factor_coul,factor_lj;
  int *ilist,*jlist,*numneigh,**firstneigh;
	int bin_intersect;
	double x_intersect;
	
  double **x = atom->x;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double *special_coul = force->special_coul;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;
	double maxcutoff,maxvalue=0.e0;
	int imax,jmax,binimax;
	double iposmax,jposmax;
	int flag=0,skip_flag=0;
  double maxx = -20.e0;
	// invoke half neighbor list (will copy or build if necessary)

  neighbor->build_one(list);

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
	
  // loop over neighbors of my atoms
  // skip if I,J are not in 2 groups
	//double *vector_allprocs = new double[maxsize];
	int sjnum = 0;
  double *vector_proc = new double[maxsize];
 	for (ii = 0; ii < maxsize; ii++) 
		vector_proc[ii] = 0.e0; 
  int count = 0;
   //std::cout<<"Inum is: "<<inum<<'\n';	
	for (ii = 0; ii < inum; ii++) {
  
	  i = ilist[ii];
	  xtmp_i = x[i][0];
    ytmp_i = x[i][1];
    ztmp_i = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];
		sjnum += jnum;
		maxcutoff = 0.e0;
  	flag = 0;	
	  for (kk = 0; kk <= atom->ntypes; kk++) {
      if(cutsq[itype][kk]> maxcutoff)
				maxcutoff=cutsq[itype][kk];
    }    
    if((ztmp_i-zpoint)*(ztmp_i-zpoint)>maxcutoff)
		{ 
				continue;
		}

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;
      jtype = type[j];
      xtmp_j = x[j][0];
      ytmp_j = x[j][1];
      ztmp_j = x[j][2];
      maxcutoff = 0.e0;
			if(xtmp_j >= maxx)
			{
				maxx = xtmp_j;
				//flag = 1;
			}
			for (kk = 0; kk <= atom->ntypes; kk++) {
      	if(cutsq[jtype][kk]> maxcutoff)
					maxcutoff=cutsq[itype][kk];
    	}   

	
    	if((ztmp_j-zpoint)*(ztmp_j-zpoint)>maxcutoff)
			{ 
				count+=1;
				continue;
			}
			 
      if (((ztmp_j >= zpoint) && (ztmp_i >= zpoint))||
					((ztmp_j < zpoint) && (ztmp_i < zpoint))) 
				continue;
				
      delx = xtmp_i - xtmp_j;
      dely = ytmp_i - ytmp_j;
      delz = ztmp_i - ztmp_j;
      rsq = delx*delx + dely*dely + delz*delz;
      
				
      if (rsq < cutsq[itype][jtype]) {
       					
				x_intersect = (zpoint - ztmp_j)*(xtmp_i - xtmp_j)/(ztmp_i - ztmp_j) + xtmp_j;
       	//Periodic BCs
				x_intersect = (x_intersect > 0) ? x_intersect : x_intersect + domain->boxhi[0];
				x_intersect = (x_intersect > domain->boxhi[0]) ? x_intersect - domain->boxhi[0] : x_intersect;
				
				bin_intersect = ceil(x_intersect/xbin_size);
				eng = pair->single(i,j,itype,jtype,rsq,factor_coul,factor_lj,fpair);
				if((bin_intersect <= maxsize)&&(bin_intersect>=1))
				{
						sgn_i = sgn_j = -1;
						if((ztmp_i - zpoint) > 0.e0)
								sgn_i = 1;
						if((ztmp_j - zpoint) > 0.e0)	
								sgn_j = 1;

				/*		if(fpair*delz*(sgn_i-sgn_j) <= maxvalue)
						{
							maxvalue = fpair*delz*(sgn_i-sgn_j);
							imax = itype;
							jmax = jtype;
							iposmax = x[i][0];
							jposmax = x[j][0];
							binimax = bin_intersect;
						}	*/
						vector_proc[bin_intersect-1] += fpair*delz*(sgn_i
					-sgn_j);
						
				}   
      }
				
			
    }
	//	if(maxx <=-15.e0)
		//	std::cout<<"Something's off! "<<jnum<<"\n";
  }
		
	//std::cout<<"Maximum x seen"<<maxx<<'\n';
	//if(maxx==1000.e0) 
		//std::cout<<"The jth atom is at position: "<<xtmp_j;
//	std::cout<<"Sum of jnum is: "<<sjnum<<'\n';	
  if(count==inum)
		std::cout<<"Something's wrong"<<inum<<'\n';
  //std::cout<<vector_proc[1]<<"here. \n";
	double vector_allprocs[maxsize];
	for(kk=0;kk<maxsize;kk++)
		vector_allprocs[kk] = 0.e0;
	MPI_Allreduce(vector_proc,vector_allprocs,maxsize,MPI_DOUBLE,MPI_SUM,world);
  
		//scalar += 0.e0;
	//std::cout<<"Itype: "<<imax<<" Jtype : "<<jmax<<" \n";
	//std::cout<<"Ith Pos: "<<iposmax<<" Jth Pos : "<<jposmax<<" \n";
	//std::cout<<"Bin intersect:"<<binimax<<"Max Value: "<<maxvalue<<"\n";
//	maxvalue = 1.0;

	int id;
	int ierr = MPI_Comm_rank(MPI_COMM_WORLD, &id); 
		
	for(kk=0; kk<maxsize; kk++)
	{
		vector[kk] = 0.5*vector_allprocs[kk];
			
	}
			
	//std::cout<<"bin index : "<<binimax<<"\n";
	//Write vector to a file.
	/*if(id==0)
	{
		std::ofstream myfile;
		myfile.open("check_stress.dat",std::ofstream::app);
		for(kk=0;kk<maxsize;kk++) {
			myfile << vector[kk] ;
			myfile << "\n";
		}
		myfile.close(); 
	}
	
	*/
	free(vector_proc);
	//Add kinetic term
	


}


