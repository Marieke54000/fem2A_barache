#pragma once

#include "mesh.h"
#include "fem.h"
#include <math.h>
#include <cmath>
#include <iostream>



namespace FEM2A {
    namespace Simu {

        //#################################
        //  Useful functions
        //#################################

        double unit_fct( vertex v )
        {
            return 1.;
        }

        double zero_fct( vertex v )
        {
            return 0.;
        }

        double xy_fct( vertex v )
        {
            return v.x + v.y;
        }
        
        double sinus_bump (vertex v )
        {
        	const double pi = std::acos(-1);
        	return 2 * (pi*pi) * std::sin(pi*v.x) * std::sin(pi*v.y);
        }

        //#################################
        //  Simulations
        //#################################
        /*
         bool test_assemble_elementary_matrix ( ){
	    // teste avec la fonction constante=1
	    
	     	    
            Mesh mesh;
       	    mesh.load("data/square.mesh");
       	    ElementMapping elt_mapping = ElementMapping ( mesh, false, 4 );
       	    ShapeFunctions ref_functions = ShapeFunctions(2,1);
       	    Quadrature quad = Quadrature::get_quadrature(2);
       	    DenseMatrix Ke;
       	    assemble_elementary_matrix(elt_mapping,ref_functions,quad,unit_fct,Ke);
       	           	    
            //Ke.print();
            
            SparseMatrix K(mesh.nb_vertices());
            local_to_global_matrix (mesh, 4, Ke, K);
            
            //K.print();
           
            
            return true;
            }*/
/*
        void pure_dirichlet_pb( const std::string& mesh_filename, bool verbose )
        {
        
            std::cout << "Solving a pure Dirichlet problem" << std::endl;
            if ( verbose ) {
                std::cout << " with lots of printed details..." << std::endl;
            }
          
            Mesh mesh;
            mesh.load (mesh_filename);
            
            // créer attribute_is_dirichlet avec Mesh::set_attribute
            std::vector< bool > attribute_is_dirichlet(2,false);
            attribute_is_dirichlet[1]=true; 
            mesh.set_attribute (unit_fct,1,true);    
            
             // créer values qui contient valeurs de g 
       	    std::vector< double > values(mesh.nb_vertices());
       	    for (int i = 0 ; i < mesh.nb_vertices() ; ++i ){
       	    	vertex v = mesh.get_vertex(i);
       	    	double g = xy_fct(v);
       	    	values[i]= g ;
       	    	}
                       
       	    //K
       	    SparseMatrix K(mesh.nb_vertices());       	    
       	    for ( int triangle = 0 ; triangle < mesh.nb_triangles() ; ++triangle){
       	    	// création de Ke associé à chaque triangle
       	    	ElementMapping elt_mapping ( mesh, false, triangle );
       	    	ShapeFunctions ref_functions (2,1);
       	    	Quadrature quad = Quadrature::get_quadrature(2);       	    
       	    	DenseMatrix Ke;
       	    	Ke.set_size (3,3);
       	    	assemble_elementary_matrix(elt_mapping,ref_functions,quad,unit_fct,Ke);
       	    	local_to_global_matrix (mesh,triangle, Ke, K);
       	    	}
       	          	    	
       	    //F
       	    std::vector< double > F(mesh.nb_vertices(),0);
       	    	
       	    // appliquer Dirichlet
       	    apply_dirichlet_boundary_conditions ( mesh, attribute_is_dirichlet, values,K,F);
       	    std :: vector <double> u(mesh.nb_vertices());
       	    solve(K,F,u);
       	    //std :: string export_name = "pureDirichletsquare";
       	    std :: string export_name = "pureDirichletsquarefine";
       	    mesh.save (export_name + ".mesh");
       	    save_solution(u, export_name+ ".bb");	
       	            	           	          	           	                
        
	}
	
    
    
     void dirichlet_pb_ts( const std::string& mesh_filename, bool verbose )
        {
        
            std::cout << "Solving a Dirichlet problem with source term" << std::endl;
            if ( verbose ) {
                std::cout << " with lots of printed details..." << std::endl;
            }
          
            Mesh mesh;
            mesh.load (mesh_filename);
            
            // créer attribute_is_dirichlet avec Mesh::set_attribute
            std::vector< bool > attribute_is_dirichlet(2,false);
            attribute_is_dirichlet[1]=true; 
            mesh.set_attribute (unit_fct,1,true);    
            
             // créer values qui contient valeurs de g 
       	    std::vector< double > values(mesh.nb_vertices(),0);
                       
       	    //K et F
       	    SparseMatrix K(mesh.nb_vertices());
       	    std::vector< double > F(mesh.nb_vertices());       	    
       	    for ( int triangle = 0 ; triangle < mesh.nb_triangles() ; ++triangle){
       	    	// création de Ke associé à chaque triangle
       	    	ElementMapping elt_mapping ( mesh, false, triangle );
       	    	ShapeFunctions ref_functions (2,1);
       	    	Quadrature quad = Quadrature::get_quadrature(2);       	    
       	    	DenseMatrix Ke ;       	    	
       	    	Ke.set_size (3,3);
       	    	assemble_elementary_matrix(elt_mapping,ref_functions,quad,unit_fct,Ke);
       	    	local_to_global_matrix (mesh,triangle, Ke, K);
       	    	
       	    	// création de Fe associé à chaque triangle
       	    	std::vector< double >Fe ;
       	    	assemble_elementary_vector(elt_mapping,ref_functions,quad,unit_fct,Fe);
       	    	local_to_global_vector(mesh,false,triangle,Fe,F);
       	    	
       	    	}
       	   
       	    // appliquer Dirichlet
       	    apply_dirichlet_boundary_conditions ( mesh, attribute_is_dirichlet, values,K,F);
       	    std :: vector <double> u(mesh.nb_vertices());
       	    solve(K,F,u);
       	    std :: string export_name = "Dirichlet_ts_square";
       	    //std :: string export_name = "Dirichlet_ts_square_fine";
       	    mesh.save (export_name + ".mesh");
       	    save_solution(u, export_name+".bb");	
       	            	           	          	           	                
        
	}
	*/
	
	/*void pb_sinus_bump( const std::string& mesh_filename, bool verbose )
        {
        
            std::cout << "Solving a Dirichlet problem with source term" << std::endl;
            if ( verbose ) {
                std::cout << " with lots of printed details..." << std::endl;
            }
          
            Mesh mesh;
            mesh.load (mesh_filename);
            
            // créer attribute_is_dirichlet avec Mesh::set_attribute
            std::vector< bool > attribute_is_dirichlet(2,false);
            attribute_is_dirichlet[1]=true; 
            mesh.set_attribute (unit_fct,1,true);    
            
             // créer values qui contient valeurs de g 
       	    std::vector< double > values(mesh.nb_vertices(),0);
                       
                                   
       	    //K et F
       	    SparseMatrix K(mesh.nb_vertices());
       	    std::vector< double > F(mesh.nb_vertices());       	    
       	    for ( int triangle = 0 ; triangle < mesh.nb_triangles() ; ++triangle){
       	    	// création de Ke associé à chaque triangle
       	    	ElementMapping elt_mapping ( mesh, false, triangle );
       	    	ShapeFunctions ref_functions (2,1);
       	    	Quadrature quad = Quadrature::get_quadrature(2);       	    
       	    	DenseMatrix Ke ;       	    	
       	    	Ke.set_size (3,3);
       	    	assemble_elementary_matrix(elt_mapping,ref_functions,quad,unit_fct,Ke);
       	    	local_to_global_matrix (mesh,triangle, Ke, K);
       	    	
       	    	// création de Fe associé à chaque triangle
       	    	std::vector< double >Fe ;
       	    	assemble_elementary_vector(elt_mapping,ref_functions,quad,sinus_bump,Fe);
       	    	local_to_global_vector(mesh,false,triangle,Fe,F);
       	    	
       	    	}
       	   
       	    // appliquer Dirichlet
       	    apply_dirichlet_boundary_conditions ( mesh, attribute_is_dirichlet, values,K,F);
       	    std :: vector <double> u(mesh.nb_vertices());
       	    solve(K,F,u);
       	    //std :: string export_name = "sinus_bump_square";
       	    std :: string export_name = "test_sinus_bump_square_fine";
       	    mesh.save (export_name + ".mesh");
       	    save_solution(u, export_name+".bb");
       	           	          		       	            	           	          	           	              
	}*/
	
	/*void sol_exacte(const std::string& mesh_filename,bool verbose){
	    Mesh mesh;
            mesh.load(mesh_filename);
            
            const double pi = std::acos(-1);
            std::vector<double> u(mesh.nb_vertices(),0);
            for (int i =0 ; i < mesh.nb_vertices() ; ++ i ){
            	u[i] = std :: sin(pi*mesh.get_vertex(i).x)*std::sin(pi*mesh.get_vertex(i).y);
            	}
            std :: string export_name = "sinus_bump_ex_square_fine";
       	    mesh.save (export_name + ".mesh");
       	    save_solution(u, export_name+".bb");
       	    }	*/
	 
	 
	 void ecart_solnum_solanalyt(const std::string& mesh_filename, bool verbose){
            Mesh mesh;
            mesh.load(mesh_filename);
            
            std::vector< double > val(mesh.nb_vertices(),0);
            std::vector< double > ecart(mesh.nb_vertices(),0);
            const double pi = std::acos(-1);
            
            // créer attribute_is_dirichlet avec Mesh::set_attribute
            std::vector< bool > attribute_is_dirichlet(2,false);
            attribute_is_dirichlet[1]=true; 
            mesh.set_attribute (unit_fct,1,true);    
            
             // créer values qui contient valeurs de g 
       	    std::vector< double > values(mesh.nb_vertices(),0);
                       
            //K et F
            SparseMatrix K(mesh.nb_vertices());
            std::vector< double > F(mesh.nb_vertices());            
            for ( int triangle = 0 ; triangle < mesh.nb_triangles() ; ++triangle){
       	    	// création de Ke associé à chaque triangle
       	    	ElementMapping elt_mapping ( mesh, false, triangle );
       	    	ShapeFunctions ref_functions (2,1);
       	    	Quadrature quad = Quadrature::get_quadrature(2);       	    
       	    	DenseMatrix Ke ;       	    	
       	    	Ke.set_size (3,3);
       	    	assemble_elementary_matrix(elt_mapping,ref_functions,quad,unit_fct,Ke);
       	    	local_to_global_matrix (mesh,triangle, Ke, K);
       	    	
       	    	// création de Fe associé à chaque triangle
       	    	std::vector< double >Fe ;
       	    	assemble_elementary_vector(elt_mapping,ref_functions,quad,sinus_bump,Fe);
       	    	local_to_global_vector(mesh,false,triangle,Fe,F);
       	    	
       	    	}
       	    	
       	    // appliquer Dirichlet
       	    apply_dirichlet_boundary_conditions ( mesh, attribute_is_dirichlet, values,K,F);
       	    std :: vector <double> u(mesh.nb_vertices());
       	    solve(K,F,u);
            
           // ecart 
           for (int i =0 ; i < mesh.nb_vertices() ; ++i ){
           ecart[i] = std::abs((std :: sin(pi*mesh.get_vertex(i).x)*std::sin(pi*mesh.get_vertex(i).y))-u[i]);
           }

           std :: string export_name = "ecart_square_fine";
       	   mesh.save (export_name + ".mesh");
       	   save_solution(ecart, export_name+".bb");


	}
	}
	

}	

