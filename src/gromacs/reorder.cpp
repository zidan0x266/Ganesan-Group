#include <stdio.h>
//#include <tchar.h>
#include <iostream>
using namespace std;
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#define real             float
#define EMPTY            -1
#define max_atomtype_CG  100                        //maximal number of CG atom types 
#define max_atomtype_AA  100                        //maximal number of AA atom types
#define max_restype      100                        //maximal number of residue types
#define max_natoms_AA    1000000                     //maximal atom number for AA atoms
#define natoms_CG        200000                      //the number of total atoms 
#define total_res_num    100000                      //the number of residues
#define scalef           0.4047                        //map the Coordinates from Ideal to Real model

//read in parameters in CG2AA library files
/*res_atom_order.dat*/
int  res_typenum;                                                //the number of residue types
int  res_atom_typenum[max_restype];                              //the number of atom types in each residue
int  res_atom_order[max_restype][max_atomtype_AA];               //the order of atoms in each residue
char res_name[max_restype][6];                                   //residue names in 'res_atom_order.dat' file                                                             
char res_AAatom_name[max_restype][max_atomtype_AA][6];           //the name AA atom in each residue

typedef struct 
{
   //coarse grained atoms and coordinate
	int   totalatom_CG;								  //total CG atom number	
    int   residuenr_CG[natoms_CG+1];				  //the residue number for each CG atom
	int   atomnr_CG[natoms_CG+1];					  //the assigned number for each CG atom
	char  title_CG[100];							  //title of the CG grofile
	char  residuename_CG[natoms_CG+1][6];			  //the residue name of CG atom
	char  atomname_CG[natoms_CG+1][6];				  //the CG atom name
	real  x_CG[natoms_CG+1][3];                       //coordinate of the CG atoms
	real  vx_CG[natoms_CG+1][3];		              //velocity of the CG atoms
	real  box_x,box_y,box_z;	                      //the box of the system 
}grofile;

typedef struct
{
	real mass[natoms_CG];
	int atomtype[natoms_CG];
	int numbonds[natoms_CG];                                        //how many bonds have been formed at this atom;  starting from 1 not from 0
	int monlinks[natoms_CG][4];                                     // which atom is linked at this atom by this bond; starting from [][1]
	int residue_num[natoms_CG];                                //residue number for each atom
	int charge_gn[natoms_CG];
	real charge[natoms_CG];
}atoms;

atoms work;
grofile grofile1;

void read_CG_gro(atoms *work, grofile *grofile1)
{
	int kk, res, tmp;
	bool flag_vec;
	char temp_str[100],temp_str_sub[16];
    char j[200],k[200],m[200];
    FILE *fp;
	fp=fopen("cg_conf.gro","r");

	//test which read style is suitable
	fgets(j,200,fp);
	fgets(k,200,fp);
	fgets(m,200,fp);
	kk=strlen(m);
	if (kk<50)
	{
		flag_vec=false;
	}else
	{
		flag_vec=true;
	}
	rewind(fp);

	//read in gro file		
	if (flag_vec==true)
	{
		fgets(grofile1->title_CG,100,fp);      //title
		fgets(temp_str,10,fp);                //totalatom
	    grofile1->totalatom_CG=atoi(temp_str);
		for (int i=1;i<=grofile1->totalatom_CG;i++)
		{
			fgets(temp_str,100,fp);                   //read whole line
			
			strncpy(temp_str_sub,temp_str,15);
			temp_str_sub[15]='\0';

			sscanf(temp_str_sub,"%i%s%s",&grofile1->residuenr_CG[i],&grofile1->residuename_CG[i],&grofile1->atomname_CG[i]);
			sscanf(temp_str+15,"%5i%8f%8f%8f%8f%8f%8f",&grofile1->atomnr_CG[i],&grofile1->x_CG[i][0],&grofile1->x_CG[i][1],&grofile1->x_CG[i][2],&grofile1->vx_CG[i][0],&grofile1->vx_CG[i][1],&grofile1->vx_CG[i][2]);
       /*
		    work->x[i][0]=grofile1->x_CG[i][0];
            work->x[i][1]=grofile1->x_CG[i][1];
			work->x[i][2]=grofile1->x_CG[i][2];
			work->vx[i][0]=grofile1->vx_CG[i][0];
			work->vx[i][1]=grofile1->vx_CG[i][1];
			work->vx[i][2]=grofile1->vx_CG[i][2];
		*/	
		}
        res = 0;
        for (int i=1;i<=grofile1->totalatom_CG;i++)
        {
            if(grofile1->residuenr_CG[i] != grofile1->residuenr_CG[i-1])
                res++;
        }
        for (int i=1;i<=grofile1->totalatom_CG;i++)
        {
            if(grofile1->residuenr_CG[i] < grofile1->residuenr_CG[i-1])
            {
                tmp = grofile1->residuenr_CG[i];
                grofile1->residuenr_CG[i] = grofile1->residuenr_CG[i-1] + 1;
                if(grofile1->residuenr_CG[i+1] == tmp)
                    grofile1->residuenr_CG[i+1] = grofile1->residuenr_CG[i];
                if(grofile1->residuenr_CG[i+2] == tmp)
                    grofile1->residuenr_CG[i+2] = grofile1->residuenr_CG[i];

            }
        }
		fscanf(fp,"%f%f%f",&grofile1->box_x,&grofile1->box_y,&grofile1->box_z);
	}
	else
	{
		fgets(grofile1->title_CG,100,fp);//title
		printf("The title of the simulated system is %s.\n", grofile1->title_CG);
		fgets(temp_str,10,fp);//totalatom
		grofile1->totalatom_CG=atoi(temp_str);
		printf("The total number of coarse-grained beads is %d.\n", grofile1->totalatom_CG);
		for (int i=1;i<=grofile1->totalatom_CG;i++)
		{
		    fgets(temp_str,100,fp);//read whole line
			strncpy(temp_str_sub,temp_str,15);
			temp_str_sub[15]='\0';
			sscanf(temp_str_sub,"%5i%5s%5s",&grofile1->residuenr_CG[i],&grofile1->residuename_CG[i],&grofile1->atomname_CG[i]);
			sscanf(temp_str+15,"%5i%8f%8f%8f",&grofile1->atomnr_CG[i],&grofile1->x_CG[i][0],&grofile1->x_CG[i][1],&grofile1->x_CG[i][2]);
            //printf("%d %f %f %f\n",grofile1->atomnr_CG[i],grofile1->x_CG[i][0],grofile1->x_CG[i][1],grofile1->x_CG[i][2] );
		 /* work->x[i][0]=grofile1->x_CG[i][0];
            work->x[i][1]=grofile1->x_CG[i][1];
			work->x[i][2]=grofile1->x_CG[i][2];	
		 */
		}
		res = 0;
		for (int i=1;i<=grofile1->totalatom_CG;i++)
        {
            if(grofile1->residuenr_CG[i] != grofile1->residuenr_CG[i-1])
                res++;
        }
        for (int i=1;i<=grofile1->totalatom_CG;i++)
        {
            if(grofile1->residuenr_CG[i] < grofile1->residuenr_CG[i-1])
            {
                tmp = grofile1->residuenr_CG[i];
                grofile1->residuenr_CG[i] = grofile1->residuenr_CG[i-1] + 1;
                if(grofile1->residuenr_CG[i+1] == tmp)
                    grofile1->residuenr_CG[i+1] = grofile1->residuenr_CG[i];
                if(grofile1->residuenr_CG[i+2] == tmp)
                    grofile1->residuenr_CG[i+2] = grofile1->residuenr_CG[i];

            }
        }
		printf("The total number of residue is %d", res);
		fscanf(fp,"%f%f%f",&grofile1->box_x,&grofile1->box_y,&grofile1->box_z);
	}
	fclose(fp);
}

void output_CG_gro(grofile *grofile1)
{
    int i,j,k;
    int res1,CG_type,re_order_AA;
    int AA_number=0;
    FILE *fp;

    printf("The title of the atimstic system is %s.\n", grofile1->title_CG);
    printf("the total number of atoms is %d.\n", grofile1->totalatom_CG);
    fp=fopen("cg_reorder.gro","w");
    fprintf(fp,"%s",grofile1->title_CG);
    fprintf(fp,"%d\n",grofile1->totalatom_CG);
    for(i=1; i<=grofile1->totalatom_CG; i++)
    {
        fprintf(fp,"%5i%-5s%5s%5i%8.3f%8.3f%8.3f\n",grofile1->residuenr_CG[i]%1000000,
                grofile1->residuename_CG[i],grofile1->atomname_CG[i],grofile1->atomnr_CG[i]%1000000,
                grofile1->x_CG[i][0]*scalef,grofile1->x_CG[i][1]*scalef,grofile1->x_CG[i][2]*scalef);
    }
    fprintf(fp,"%8.5f %8.5f %8.5f\n",grofile1->box_x*scalef,grofile1->box_y*scalef,grofile1->box_z*scalef);
    fclose(fp);
    printf("the cg_reorder.gro is finished.\n");
}

int main(int argc,char *argv[])
{
	read_CG_gro(&work, &grofile1);
    output_CG_gro(&grofile1);
	return 0;
}
