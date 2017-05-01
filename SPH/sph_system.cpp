/** File:		sph_system.cpp
** Author:		Dongli Zhang
** Contact:	dongli.zhang0129@gmail.com
**
** Copyright (C) Dongli Zhang 2013
**
** This program is free software;  you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation; either version 2 of the License, or
** (at your option) any later version.
**
** This program is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY;  without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See
** the GNU General Public License for more details.
**
** You should have received a copy of the GNU General Public License
** along with this program;  if not, write to the Free Software
** Foundation, 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
*/

#include "sph_system.h"
#include "sph_header.h"

SPHSystem::SPHSystem()
{
	numIterations = 5;
	max_particle = 30000;
	num_particle = 0;

	kernel = 0.1f;//0.1f
	mass =1.0f;

	world_size.x = 0.64f;
	world_size.y = 0.64f;
	world_size.z = 0.64f;
	cell_size = kernel;
	grid_size.x = (uint)ceil(world_size.x / cell_size);
	grid_size.y = (uint)ceil(world_size.y / cell_size);
	grid_size.z = (uint)ceil(world_size.z / cell_size);
	tot_cell = grid_size.x*grid_size.y*grid_size.z;

	gravity.x = 0.0f;
	gravity.y = -9.8f;
	gravity.z = 0.0f;
	wall_damping = -0.5f;
	rest_density = 6000.0f;


	time_step = 0.003f;//Ê±¼ä²½³¤£¿


	poly6_value = 315.0f / (64.0f * PI * pow(kernel, 9));;
	spiky_value = 45.0f / (PI * pow(kernel, 6));
	/*visco_value = 45.0f / (PI * pow(kernel, 6));*/

	K = 0.00001f;
	C = 0.0025f; //0.01f;
	dqMag = 0.2f * kernel;
	wQH = poly6_value * pow((kernel *kernel - dqMag * dqMag), 3);



	kernel_2 = kernel*kernel;
	self_dens = mass*poly6_value*pow(kernel, 6);//×ÔÉíÃÜ¶È

	lambdaEps = 600.0f;//

	mem = (Particle *)malloc(sizeof(Particle)*max_particle);
	cell = (Particle **)malloc(sizeof(Particle *)*tot_cell);
	

	sys_running = 0;

	printf("Initialize SPH:\n");
	printf("World Width : %f\n", world_size.x);
	printf("World Height: %f\n", world_size.y);
	printf("World Length: %f\n", world_size.z);
	printf("Cell Size  : %f\n", cell_size);
	printf("Grid Width : %u\n", grid_size.x);
	printf("Grid Height: %u\n", grid_size.y);
	printf("Grid Length: %u\n", grid_size.z);
	printf("Total Cell : %u\n", tot_cell);
	printf("Poly6 Kernel: %f\n", poly6_value);
	printf("Spiky Kernel: %f\n", spiky_value);

	printf("Self Density: %f\n", self_dens);
}

SPHSystem::~SPHSystem()
{
	free(mem);
	free(cell);
}

inline void operator +=(float3 &a, float3 b)
{
	a.x += b.x;
	a.y += b.y;
	a.z += b.z;
}

inline  float length(float3 v)
{
	return sqrtf(v.x * v.x + v.y * v.y + v.z * v.z);
}

inline float3 make_float3(float x, float y, float z)
{
	float3 t; t.x = x; t.y = y; t.z = z; return t;
}

inline float3 operator*(float3 a, float b)
{
	return make_float3(a.x * b, a.y * b, a.z * b);
}

inline float3 operator-(float3 a, float3 b)
{
	return make_float3(a.x - b.x, a.y - b.y, a.z - b.z);
}

inline float3 operator/(float3 a, float b)
{
	return make_float3(a.x / b, a.y / b, a.z / b);
}
inline void operator*=(float3 &a, float b)
{
	a.x *= b;
	a.y *= b;
	a.z *= b;
}




float SPHSystem::sCorrCalc(Particle *p, Particle *np)
{
	//Get Density from WPoly6
	float3 r;
	r.x = p->pos.x - np->pos.x;
	r.y = p->pos.y - np->pos.y;
	r.z = p->pos.z - np->pos.z;
	float len = length(r);
	float corr = poly6_value*pow((kernel_2 - pow(len, 2)), 3) / wQH;
	corr *= corr * corr * corr;
	return -K * corr;//
}




void SPHSystem::animation()
{
	if (sys_running == 0)
	{
		return;
	}
	neighbour.resize(max_particle);//
    build_table();
	for (int i = 0; i < numIterations; i++)
	{
		predictPositions();
	
			calcLambda();

			calcDeltaP();
	
			applyDeltaP();
		
	}
	
	updateVelocities();
	/*printf("do");*/
	updateXSPHVelocities();
	/*printf("do");*/
	neighbour.clear();
}

void SPHSystem::init_system()
{
	float3 pos;
	float3 vel;

	vel.x = 0.0f;
	vel.y = 0.0f;
	vel.z = 0.0f;

	for (pos.x = world_size.x*0.0f; pos.x<world_size.x*0.6f; pos.x += (kernel*0.5f))
	{
		for (pos.y = world_size.y*0.0f; pos.y<world_size.y*0.9f; pos.y += (kernel*0.5f))
		{
			for (pos.z = world_size.z*0.0f; pos.z<world_size.z*0.6f; pos.z += (kernel*0.5f))
			{
				add_particle(pos, vel);
			}
		}
	}

	printf("Init Particle: %u\n", num_particle);
}

void SPHSystem::add_particle(float3 pos, float3 vel)
{
	Particle *p = &(mem[num_particle]);

	p->id = num_particle;
	
	p->deltaPs=make_float3(0.0f ,0.0f ,0.0f);
	
	
	p->lambda= 0.0f;

	p->oldPos = pos;
	p->pos = pos;
	p->vel = vel;


	p->dens = 0.0f;
    p->visc = make_float3(0.0f, 0.0f, 0.0f);

	p->next = NULL;

	num_particle++;
}

void SPHSystem::build_table()//½¨Á¢Ë÷Òý
{
	Particle *p;
	uint hash;

	for (uint i = 0; i<tot_cell; i++)
	{
		cell[i] = NULL;
	}

	for (uint i = 0; i<num_particle; i++)
	{
		p = &(mem[i]);
		hash = calc_cell_hash(calc_cell_pos(p->pos));

		if (cell[hash] == NULL)
		{
			p->next = NULL;
			cell[hash] = p;
		}
		else
		{
			p->next = cell[hash];
			cell[hash] = p;
		}
	}
}

void SPHSystem::comp_dens()
{
	Particle *p;
	Particle *np;

	int3 cell_pos;
	int3 near_pos;
	uint hash;

	float3 rel_pos;
	float r2;

	for (uint i = 0; i<num_particle; i++)
	{
		p = &(mem[i]);
		cell_pos = calc_cell_pos(p->pos);

	    //ÁÚÓòËÑË÷
		for (int x = -1; x <= 1; x++)   
		{
			for (int y = -1; y <= 1; y++)
			{
				for (int z = -1; z <= 1; z++)
				{
					near_pos.x = cell_pos.x + x;
					near_pos.y = cell_pos.y + y;
					near_pos.z = cell_pos.z + z;
					hash = calc_cell_hash(near_pos);

					if (hash == 0xffffffff)
					{
						continue;
					}

					np = cell[hash];
					while (np != NULL)
					{
						rel_pos.x = np->pos.x - p->pos.x;
						rel_pos.y = np->pos.y - p->pos.y;
						rel_pos.z = np->pos.z - p->pos.z;
						r2 = rel_pos.x*rel_pos.x + rel_pos.y*rel_pos.y + rel_pos.z*rel_pos.z;

						if (r2<INF || r2 >= kernel_2)
						{
							np = np->next;
							continue;
						}

						p->dens = p->dens + mass * poly6_value * pow(kernel_2 - r2, 3);

   						neighbour[i].push_back(*np);
						/*printf("do");*/
						np = np->next;
					}
				}
			}
		}

		//p->dens = p->dens + self_dens;//pbd中没加self-dens

	}
}

void SPHSystem::calcLambda()
{
	Particle *p;
	

	int3 cell_pos;
	int3 near_pos;
	uint hash;

	float3 rel_pos;
	float r2;

	for (uint i = 0; i < num_particle; i++)
	{
		p = &(mem[i]);
		
		float densityConstraint = (p->dens / rest_density) - 1;//pbd restdesity=6378

		float3 gradientI = make_float3(0.0f, 0.0f, 0.0f);
		float sumGradients = 0.0f;
		
					for(uint j = 0;j < neighbour[i].size();j++)
					{
						rel_pos.x = neighbour[i][j].pos.x - p->pos.x;
						rel_pos.y = neighbour[i][j].pos.y - p->pos.y;
						rel_pos.z = neighbour[i][j].pos.z - p->pos.z;
						r2 = rel_pos.x*rel_pos.x + rel_pos.y*rel_pos.y + rel_pos.z*rel_pos.z;

						
						float r = sqrt(r2);
						float coeff = (kernel - r)*(kernel - r);
						coeff *= spiky_value;
						coeff /=-r;
						float3 gradientJ = rel_pos*coeff;
						gradientI += gradientJ;
										
					}
		
		//Add the particle i gradient magnitude squared to sum
		sumGradients += pow(length(gradientI), 2);
		
		p->lambda = (-1 * densityConstraint) / (sumGradients + lambdaEps);
		/*printf("%f/n", p->lambda);*/
	}
}

void SPHSystem::predictPositions()
{
	Particle *p;
	for (uint i = 0; i < num_particle; i++)
	{
		p = &(mem[i]);

		p->vel.x = p->vel.x + gravity.x*time_step;
		p->vel.y = p->vel.y + gravity.y*time_step;
		p->vel.z = p->vel.z + gravity.z*time_step;

		p->pos.x = p->pos.x + p->vel.x *time_step;
		p->pos.y = p->pos.y + p->vel.y*time_step;
		p->pos.z = p->pos.z + p->vel.z*time_step;
		/*printf("%f/n", p->pos.y);*/
		if (p->pos.x >= world_size.x - BOUNDARY)
		{
			p->vel.x = p->vel.x*wall_damping;
			p->pos.x = world_size.x - BOUNDARY;
		}

		if (p->pos.x < 0.0f)
		{
			p->vel.x = p->vel.x*wall_damping;
			p->pos.x = 0.0f;
		}

		if (p->pos.y >= world_size.y - BOUNDARY)
		{
			p->vel.y = p->vel.y*wall_damping;
			p->pos.y = world_size.y - BOUNDARY;
		}

		if (p->pos.y < 0.0f)
		{
			p->vel.y = p->vel.y*wall_damping;
			p->pos.y = 0.0f;
		}

		if (p->pos.z >= world_size.z - BOUNDARY)
		{
			p->vel.z = p->vel.z*wall_damping;
			p->pos.z = world_size.z - BOUNDARY;
		}

		if (p->pos.z < 0.0f)
		{
			p->vel.z = p->vel.z*wall_damping;
			p->pos.z = 0.0f;
		}
	}
	
}

void SPHSystem::calcDeltaP()
{
	Particle *p;

	float3 rel_pos;
	float r2;



	for (uint i = 0; i < num_particle; i++)
	{
		p = &(mem[i]);

		float densityConstraint = (p->dens / rest_density) - 1;//pbd restdesity=6378
		p->deltaPs = make_float3(0.0f, 0.0f, 0.0f);
		float3 deltaP = make_float3(0.0f, 0.0f, 0.0f);
		

		for (uint j = 0; j<neighbour[i].size(); j++)
		{
			rel_pos.x = neighbour[i][j].pos.x - p->pos.x;
			rel_pos.y = neighbour[i][j].pos.y - p->pos.y;
			rel_pos.z = neighbour[i][j].pos.z - p->pos.z;
			r2 = rel_pos.x*rel_pos.x + rel_pos.y*rel_pos.y + rel_pos.z*rel_pos.z;


			float r = sqrt(r2);

			float lambdaSum = p->lambda + neighbour[i][j].lambda;
						
			float sCorr = sCorrCalc(p, &neighbour[i][j]);
			float coeff = (kernel - r)*(kernel - r);
			coeff *= spiky_value;
			coeff /= -r;
			deltaP += rel_pos*coeff* (lambdaSum + sCorr);
						
						
		}
			
		p->deltaPs = deltaP /rest_density;
		
	}

}


void SPHSystem::applyDeltaP()
{
	Particle *p;
	for (uint i = 0; i < num_particle; i++)
	{
		p = &(mem[i]);
		p->oldPos = p->pos;
		p->pos += p->deltaPs;
		if (p->pos.x >= world_size.x - BOUNDARY)
		{
			p->vel.x = p->vel.x*wall_damping;
			p->pos.x = world_size.x - BOUNDARY;
		}

		if (p->pos.x < 0.0f)
		{
			p->vel.x = p->vel.x*wall_damping;
			p->pos.x = 0.0f;
		}

		if (p->pos.y >= world_size.y - BOUNDARY)
		{
			p->vel.y = p->vel.y*wall_damping;
			p->pos.y = world_size.y - BOUNDARY;
		}

		if (p->pos.y < 0.0f)
		{
			p->vel.y = p->vel.y*wall_damping;
			p->pos.y = 0.0f;
		}

		if (p->pos.z >= world_size.z - BOUNDARY)
		{
			p->vel.z = p->vel.z*wall_damping;
			p->pos.z = world_size.z - BOUNDARY;
		}

		if (p->pos.z < 0.0f)
		{
			p->vel.z = p->vel.z*wall_damping;
			p->pos.z = 0.0f;
		}
		
	}

}



void SPHSystem::updateVelocities()
{
	Particle *p;
	
	float3 rel_pos;
	float r2;

for (uint i = 0; i < num_particle; i++)
{
	p = &(mem[i]);
    p->vel = (p->pos - p->oldPos) / time_step;//
	float densityConstraint = (p->dens / rest_density) - 1;//pbd restdesity=6378
	p->deltaPs = make_float3(0.0f, 0.0f, 0.0f);
	float3 deltaP = make_float3(0.0f, 0.0f, 0.0f);


	for (uint j = 0; j<neighbour[i].size(); j++)
	{
		rel_pos.x = neighbour[i][j].pos.x - p->pos.x;
		rel_pos.y = neighbour[i][j].pos.y - p->pos.y;
		rel_pos.z = neighbour[i][j].pos.z - p->pos.z;
		r2 = rel_pos.x*rel_pos.x + rel_pos.y*rel_pos.y + rel_pos.z*rel_pos.z;


		float r = sqrt(r2);

		float3 velocityDiff = neighbour[i][j].vel - p->vel;
		velocityDiff *= poly6_value * pow((kernel_2 - r2), 3);
		p->visc += velocityDiff;						
						                  
	}
					
		p->deltaPs =p->visc * C;
		p->oldPos = p->pos;

	}


}


void SPHSystem::updateXSPHVelocities()
{
	Particle *p;
	for (uint i = 0; i < num_particle; i++)
	{
		p = &(mem[i]);
		p->vel += p->deltaPs* time_step;
		/*printf("%f/n", p -> vel);*/
	}
}




int3 SPHSystem::calc_cell_pos(float3 p)
{
	int3 cell_pos;
	cell_pos.x = int(floor((p.x) / cell_size));
	cell_pos.y = int(floor((p.y) / cell_size));
	cell_pos.z = int(floor((p.z) / cell_size));

	return cell_pos;
}

uint SPHSystem::calc_cell_hash(int3 cell_pos)
{
	if (cell_pos.x<0 || cell_pos.x >= (int)grid_size.x || cell_pos.y<0 || cell_pos.y >= (int)grid_size.y || cell_pos.z<0 || cell_pos.z >= (int)grid_size.z)
	{
		return (uint)0xffffffff;
	}

	cell_pos.x = cell_pos.x & (grid_size.x - 1);
	cell_pos.y = cell_pos.y & (grid_size.y - 1);
	cell_pos.z = cell_pos.z & (grid_size.z - 1);

	return ((uint)(cell_pos.z))*grid_size.y*grid_size.x + ((uint)(cell_pos.y))*grid_size.x + (uint)(cell_pos.x);
}
