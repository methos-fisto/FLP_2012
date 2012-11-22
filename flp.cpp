#include "flp.h"

std::vector<int>** pretraitement(int nbclient, int nbusine, int** mat_distance, int& nbSij){

std::vector<int>** mat_sij;
nbSij = 0;
mat_sij = new std::vector<int>*[nbclient];
for(int i = 0; i< nbclient; i++){
	mat_sij[i] = new std::vector<int>[nbusine];
}

for(int i = 0; i< nbclient; i++){
	for(int j = 0; j< nbusine; j++){
		for(int k = 0; k < nbusine; k++){
			if(mat_distance[i][k] < mat_distance[i][j]){
				mat_sij[i][j].push_back(k);
				nbSij++;
			}
		}
	}
}

return mat_sij;

}

void flp_parser(std::string fichier_arg, int**& cout_livr, int**& distance, int*& opening, int& nbClients, int& nbFacility){
	//lecture du fichier
	int i, j;
	std::ifstream* fichier_traite;
	fichier_traite = new std::ifstream(fichier_arg.c_str());
	*fichier_traite >> nbClients;
	*fichier_traite >> nbFacility;
	cout_livr = new int*[nbClients];
	distance = new int*[nbClients];
	opening = new int[nbFacility];
	for(i=0; i< nbClients; i++){
	cout_livr[i] = new int[nbFacility];
	distance[i] = new int[nbFacility];
	}
	for(i=0; i< nbClients; i++){
		for(j=0; j< nbFacility; j++){
		*fichier_traite >> cout_livr[i][j];
		}
	}
	for(i=0; i< nbClients; i++){
		for(j=0; j< nbFacility; j++){
		*fichier_traite >> distance[i][j];
		}
	}
	for(i=0; i< nbFacility; i++){
		*fichier_traite >> opening[i];
	}
	
	delete fichier_traite;
}

void flp_recursif(Result*& solopti,glp_prob* probref, int* ia , int* ja ,double* ar, int tail_mat , int nbClients, int nbFacility, int& appels){

	
    glp_smcp parm;
    glp_init_smcp(&parm);
    parm.msg_lev = GLP_MSG_OFF; /* Paramtre de GLPK dans la rsolution d'un PL en variables continues afin de couper les affichages (dans lesquels on se noierait) */
	int i;
    glp_iocp parmip;
    glp_init_iocp(&parmip);
    parmip.msg_lev = GLP_MSG_OFF;
	

        glp_load_matrix(probref,tail_mat-1,ia,ja,ar); /* chargement de la matrice dans le problme */
		//resolution
        glp_simplex(probref,&parm);
        glp_intopt(probref,&parmip);
		appels++;
		
		double* x = new double[nbFacility];
		double* y = new double[nbClients*nbFacility];
		
		int plusfaible = nbFacility+1;
		bool admissible = true;
		double z = glp_mip_obj_val(probref);
		if( z!=0 ){ //si réalisable
		if(z < solopti->val()) //si interressant
		{	
			//recuperation des valeurs
			for(i=0;i<nbFacility;i++) x[i] = glp_mip_col_val(probref, i+1);
			for(i=nbFacility;i<nbFacility + nbClients*nbFacility ;i++) y[i-nbFacility] = glp_mip_col_val(probref,i+1);
			for(i=0;i<nbFacility;i++){
				if(x[i]< 0.0001f){
					x[i]=0;
				}
			}
			for(i=0;i< nbFacility; i++){
				if((x[i] != 0 ) && (x[i] != 1)){
					admissible = false;
					if(plusfaible != nbFacility+1){
						if(x[i] < x[plusfaible]){
						plusfaible= i;
						}
					}else{
						plusfaible= i;
					}
				}
			}
			//si admissible
			if(admissible){
				//nouvelle valeur
				solopti->set_val(z);
				solopti->set_sol(x);
				solopti->set_sol2(y);
			}else{
				//sinon branch
				glp_prob *prob_temp;
				prob_temp = glp_create_prob();
				glp_copy_prob(prob_temp, probref, GLP_OFF);
				int* temp_ia = new int[tail_mat+1];
				int* temp_ja = new int[tail_mat+1];
				double* temp_ar = new double[tail_mat+1];
				for(i = 0; i< tail_mat; i++){
					temp_ia[i]=ia[i];
					temp_ja[i]=ja[i];
					temp_ar[i]=ar[i];
				}
				temp_ia[tail_mat] = glp_add_rows(prob_temp, 1);;
				temp_ja[tail_mat] = plusfaible + 1;
				temp_ar[tail_mat] = 1.0f;
				glp_set_row_bnds(prob_temp, temp_ia[tail_mat], GLP_FX, 0.0f, 0.0f);
				//branche de gauche
				flp_recursif(solopti,prob_temp, temp_ia , temp_ja , temp_ar, tail_mat + 1 , nbClients, nbFacility, appels);
				glp_set_row_bnds(prob_temp, temp_ia[tail_mat], GLP_FX, 1.0f, 1.0f);
				//branche de droite
				flp_recursif(solopti,prob_temp, temp_ia , temp_ja , temp_ar, tail_mat + 1 , nbClients, nbFacility, appels);
				delete[] x;
			delete[] y;
			delete[] temp_ia;
			delete[] temp_ja;
			delete[] temp_ar;
			glp_delete_prob(prob_temp);
			}
			
		}else{
			
			delete[] x;
			delete[] y;
			
		}}else{
			delete[] x;
			delete[] y;
		}
		
		

}

void flp_solve(const std::string fichier_arg){

	int** distance;
int** 	cout_livr;
int*  opening;
	int nbClients;
	int nbSij;
	int nbFacility;
	int i,j;
	int appels = 0;
        flp_parser(fichier_arg,cout_livr, distance, opening, nbClients ,nbFacility);  
    
	std::vector<int>** Sij = pretraitement(nbClients, nbFacility, distance, nbSij);
	
	////////////////////////
	//     code GLPK      //
	////////////////////////
	glp_prob *prob; /* dclaration du pointeur sur le problme */
    prob = glp_create_prob(); /* allocation mmoire pour le problme */
    glp_set_prob_name(prob, "FLP"); /* affectation d'un nom */
    glp_set_obj_dir(prob, GLP_MIN); /* Il s'agit d'un problme de minimisation */
    //Result* resultat;
    glp_smcp parm;
    glp_init_smcp(&parm);
    parm.msg_lev = GLP_MSG_OFF; /* Paramtre de GLPK dans la rsolution d'un PL en variables continues afin de couper les affichages (dans lesquels on se noierait) */

    glp_iocp parmip;
    glp_init_iocp(&parmip);
    parmip.msg_lev = GLP_MSG_OFF;
	
	
	glp_add_rows(prob, nbClients + (2*nbClients*nbFacility));
	for(j=1;j<=nbClients ;j++){// s x y
	glp_set_row_bnds(prob, j, GLP_FX, 1.0, 1.0);
	}
	for(j=nbClients+1;j<= (nbClients + nbClients*nbFacility);j++){// s x y
	glp_set_row_bnds(prob, j, GLP_UP, 0.0, 0.0);
	}
	for(j=(nbClients + 1 + nbClients*nbFacility);j<=nbClients + (2*nbClients*nbFacility);j++){// s x y
	glp_set_row_bnds(prob, j, GLP_UP, 0.0, 0.0);
	}
	
	glp_add_cols(prob, nbFacility+ nbClients*nbFacility);
	for(j=1;j<=nbFacility+ nbClients*nbFacility;j++){// s x y
	glp_set_col_bnds(prob, j, GLP_DB, 0.0, 1.0);
	}
	
	for(j=1;j<=nbFacility;j++){// s x y
	glp_set_obj_coef(prob, j, (double) opening[j-1] );
	}
	for(i= 1; i<= nbClients ; i++){
		for(j= 1; j<= nbFacility ; j++){
			glp_set_obj_coef(prob, nbFacility + (i-1)*nbFacility + j, (double) cout_livr[i-1][j-1]);
		}
	}
	int* ia = new int[5*nbClients*nbFacility + nbSij +1];
	int* ja = new int[5*nbClients*nbFacility + nbSij +1];
	double* ar = new double[5*nbClients*nbFacility + nbSij +1];
	
	int compteur =0 ;
	for(int i = 1; i<= nbClients; i++){//y premiere contrainte
		for(j= 1; j<= nbFacility; j++){
			compteur++;
			ia[compteur] = i;
			ja[compteur] = nbFacility + (i-1)*nbFacility + j;
			ar[compteur] = 1.0f;
		}
	}
	for(int i = 1; i<= nbClients; i++){//seconde contrainte
		for(j= 1; j<= nbFacility; j++){
			compteur++;
			ia[compteur] = nbClients + (i-1)*nbFacility + j;
			ja[compteur] = nbFacility + (i-1)*nbFacility + j;
			ar[compteur] = 1.0f;
			
			compteur++;
			ia[compteur] = nbClients + (i-1)*nbFacility + j;
			ja[compteur] = j;
			ar[compteur] = -1.0f;
		}
	}
	for(int i = 1; i<= nbClients; i++){//derniere contrainte
		for(j= 1; j<= nbFacility; j++){
		compteur++;
		ia[compteur] = nbClients + nbClients*nbFacility + (i-1)*nbFacility + j;
		ja[compteur] = j;
		ar[compteur] = 1.0f;
		
		compteur++;
		ia[compteur] = nbClients + nbClients*nbFacility + (i-1)*nbFacility + j;
		ja[compteur] = nbFacility + (i-1)*nbFacility + j;
		ar[compteur] = -1.0f;
		
		//Sij
		
		for(int k = 0; k < Sij[i-1][j-1].size(); k++){
			compteur++;
			ia[compteur] = nbClients + nbClients*nbFacility + (i-1)*nbFacility + j;
			ja[compteur] = Sij[i-1][j-1].at(k) + 1;
			ar[compteur] = -1.0f;
		}
		
		
		}
	}
		
	Result* solopti = new Result(9000000.0, new double[4], new double[4]);
	//lancement de la resolution
	flp_recursif(solopti, prob , ia , ja , ar, 5*nbClients*nbFacility + nbSij +1 , nbClients, nbFacility, appels);
	
	//affichage
	std::cout << "Nombre de résolutions : " << appels << "\n";
	std::cout <<"valeur optimale : " << solopti->val() << "\n\n" << std::endl << "Plan d'ouverture : " << "\n";
	double* x = solopti->solution();
	for(int i = 0 ; i < nbFacility; i++){
		std::cout << "x[" << i +1 << "]= " << x[i] << "\n";
	}
	std::cout <<"livraison : \n";
	double* y = solopti->solution2();
	for(int i = 0 ; i < nbClients; i++){
		for(int j = 0; j < nbFacility; j++){
		std::cout << "y[" << i +1 << "][" << j+1 << " = " << y[i*nbFacility + j] << "\n";
	}}
	
}

int main(int argc, char* argv[]){

	if(argc < 2){
		flp_solve("10-5.txt");
	}else{
		flp_solve(argv[1]);
	}

}