#include "fem.h"
#include "mesh.h"

#include <iomanip>
#include <iostream>
#include <cmath>
#include <limits>
#include <stdlib.h>
#include <assert.h>

namespace FEM2A {

    void print( const std::vector<double>& x )
    {
        for ( int i = 0; i < x.size(); ++i ) {
            std::cout << x[i] << " ";
        }
        std::cout << std::endl;
    }

    /****************************************************************/
    /* Implementation of Quadrature */
    /****************************************************************/
    int Quadrature::nb_points() const
    {
        return wxy_.size() / 3 ;
    }

    vertex Quadrature::point( int i ) const
    {
        assert( i < nb_points() ) ;
        vertex v ;
        v.x = wxy_[3 * i + 1] ;
        v.y = wxy_[3 * i + 2] ;
        return v ;
    }

    double Quadrature::weight( int i ) const
    {
        assert( i < nb_points() ) ;
        return wxy_[3 * i + 0] ;
    }

    const double triangle_P0[3] = {
        0.5, 0.333333333333333, 0.333333333333333
    };

    const double triangle_P2[9] = {
        0.166666666666667, 0.166666666666667, 0.166666666666667,
        0.166666666666667, 0.166666666666667, 0.666666666666667,
        0.166666666666667, 0.666666666666667, 0.166666666666667
    };

    const double triangle_P4[18] = {
        0.0549758718276609, 0.0915762135097707, 0.0915762135097707,
        0.0549758718276609, 0.0915762135097707, 0.816847572980459,
        0.0549758718276609, 0.816847572980459, 0.0915762135097707,
        0.111690794839006, 0.445948490915965, 0.445948490915965,
        0.111690794839006, 0.445948490915965, 0.10810301816807,
        0.111690794839006, 0.10810301816807, 0.445948490915965
    };

    const double triangle_P6[36] = {
        0.0254224531851034, 0.0630890144915022, 0.0630890144915022,
        0.0254224531851034, 0.0630890144915022, 0.873821971016996,
        0.0254224531851034, 0.873821971016996, 0.0630890144915022,
        0.0583931378631897, 0.24928674517091, 0.24928674517091,
        0.0583931378631897, 0.24928674517091, 0.501426509658179,
        0.0583931378631897, 0.501426509658179, 0.24928674517091,
        0.0414255378091868, 0.0531450498448169, 0.310352451033784,
        0.0414255378091868, 0.310352451033784, 0.0531450498448169,
        0.0414255378091868, 0.0531450498448169, 0.636502499121399,
        0.0414255378091868, 0.636502499121399, 0.0531450498448169,
        0.0414255378091868, 0.310352451033784, 0.636502499121399,
        0.0414255378091868, 0.636502499121399, 0.310352451033784
    };

    const double segment_P0[2] = {
        1., 0.5
    };

    const double segment_P2[4] = {
        0.5, 0.21132486540518708,
        0.5, 0.7886751345948129
    };

    Quadrature Quadrature::get_quadrature( int order, bool border )
    {
        double* pts = NULL;
        int nb_pts = 0;
        Quadrature Q;
        if ( order == 0 && !border ) {
            pts = const_cast<double*>(triangle_P0);
            nb_pts = 1;
        } else if ( order == 2 && !border ) {
            pts = const_cast<double*>(triangle_P2);
            nb_pts = 3;
        } else if ( order == 4 && !border ) {
            pts = const_cast<double*>(triangle_P4);
            nb_pts = 6;
        } else if ( order == 6 && !border ) {
            pts = const_cast<double*>(triangle_P6);
            nb_pts = 12;
        } else if ( order == 0 && border ) {
            pts = const_cast<double*>(segment_P0);
            nb_pts = 1;
        } else if ( order == 2 && border ) {
            pts = const_cast<double*>(segment_P2);
            nb_pts = 2;
        } else {
            std::cout << "Quadrature not implemented for order " << order << std::endl;
            assert( false );
        }
        Q.wxy_.resize(nb_pts * 3);
        for ( int i = 0; i < nb_pts; ++i ) {
            if ( !border ) {
                Q.wxy_[3*i+0] = pts[3*i+0];
                Q.wxy_[3*i+1] = pts[3*i+1];
                Q.wxy_[3*i+2] = pts[3*i+2];
            } else {
                Q.wxy_[3*i+0] = pts[2*i+0];
                Q.wxy_[3*i+1] = pts[2*i+1];
                Q.wxy_[3*i+2] = 0.;
            }
        }
        return Q;
    }

    /****************************************************************/
    /* Implementation of ElementMapping */
    /****************************************************************/
    ElementMapping::ElementMapping( const Mesh& M, bool border, int i ) // méthode constructeur 
        : border_( border )
    {
        //std::cout << "[ElementMapping] constructor for element " << i << " ";
        //std :: cout << '\n';
        if ( border ){
        	for ( int v = 0 ; v < 2 ; ++v ) vertices_.push_back (M.get_edge_vertex (i,v));
        	//for ( int v = 0 ; v < 2 ; ++v ) std :: cout << vertices_[v].x << " " << vertices_[v].y << std :: endl ;
        }
        else {
        	for ( int v = 0 ; v < 3 ; ++v ) vertices_.push_back(M.get_triangle_vertex(i,v));
        	//for ( int v = 0 ; v < 3; ++v ) std :: cout << vertices_[v].x << " " << vertices_[v].y << std :: endl ;
        }
        
    }

    vertex ElementMapping::transform( vertex x_r ) const
    {
        //std::cout << "[ElementMapping] transform reference to world space" << '\n';
        vertex r ;
        if ( border_ ){
        	double phi0 = 1 - x_r.x ;
        	double phi1 = x_r.x ;
        	r.x = phi0 * vertices_[0].x+ phi1 * vertices_[1].x;
        	r.y = phi0 * vertices_[0].y+ phi1 * vertices_[1].y;
       }
        else {
        	double phi0 = 1 - x_r.x - x_r.y ;
        	double phi1 = x_r.x ;
        	double phi2 = x_r.y ; 
        	r.x = phi0 * vertices_[0].x + phi1 * vertices_[1].x + phi2 * vertices_[2].x;
        	r.y = phi0 * vertices_[0].y + phi1 * vertices_[1].y + phi2 * vertices_[2].y;
        	
       }
        return r ;
    }

    DenseMatrix ElementMapping::jacobian_matrix( vertex x_r ) const
    {
        //std::cout << "[ElementMapping] compute jacobian matrix" << '\n';
        DenseMatrix J ;
         if ( border_ ){
         	J.set_size (2 , 1);     	
        	J.set( 0 , 0 , - vertices_[0].x + vertices_[1].x);
        	J.set( 0 , 1 , - vertices_[0].y + vertices_[1].y);
       }
        else {
        	J.set_size (2 , 2);        	
        	J.set( 0 , 0 , - vertices_[0].x +  vertices_[1].x);
        	J.set( 0 , 1 , - vertices_[0].x +  vertices_[2].x);
        	J.set ( 1 , 0 , - vertices_[0].y +  vertices_[1].y);
        	J.set ( 1 , 1 , - vertices_[0].y +  vertices_[2].y);
       }
        return J ;
    }

    double ElementMapping::jacobian( vertex x_r ) const
    {
        //std::cout << "[ElementMapping] compute jacobian determinant" << '\n';
        DenseMatrix J = jacobian_matrix (x_r);
        double determinant;
        if ( border_ ){        	
         	determinant = sqrt ( J.get (0,0) * J.get (0,0) + J.get (0,1) * J.get (1,0));
       }
        else {        	
        	determinant = J.det_2x2();
       }
       
        return determinant;
        
    }

    /****************************************************************/
    /* Implementation of ShapeFunctions */
    /****************************************************************/
    ShapeFunctions::ShapeFunctions( int dim, int order )
        : dim_( dim ), order_( order )
    {
        //std::cout << "[ShapeFunctions] constructor in dimension " << dim << '\n';
    }
        /*
        bool SF_construct = true ;
        if ( dim! = 1 && dim! = 2){
        	std::cout << "ShapeFunctions are only implemented in 1D or 2D. " << std::endl;
        	SF_construct=false;
        }
        if (order != 1 ){
        	std::cout<< "Only order-1 ShapeFunctions are implemented " << std ::endl;
        	SF_construct = false;
        }
        assert (SF_construct);
    }*/

    int ShapeFunctions::nb_functions() const
    {
        //std::cout << "[ShapeFunctions] number of functions" << '\n';
        int nb_functions ;
        if ( dim_ == 1 ){
        	nb_functions = 2 ;
        }
        if( dim_ == 2){
        	nb_functions = 3 ;
        }
       
        return nb_functions ;
    }

    double ShapeFunctions::evaluate( int i, vertex x_r ) const
  {
        //std::cout << "[ShapeFunctions] evaluate shape function " << i << '\n';
        double xi = x_r.x;
        double eta = x_r.y;
        if( dim_ == 1){
        	switch (i){
        		case (0):
        			return(1- xi);
   
        		case (1) :
        			return (xi);
        			
        			}
        }
        if ( dim_ == 2 ){
        	switch (i){
        		case (0) :
        			return (1- xi-eta);
        			
        		case (1) :
        			return (xi);
        			
        		case (2) :
        			return (eta);
        			
        			}
        } 
        return 0.;       			        	        
    }

    vec2 ShapeFunctions::evaluate_grad( int i, vertex x_r ) const
    {
        //std::cout << "[ShapeFunctions] evaluate gradient shape function " << i << '\n';
        vec2 g ;
        if( dim_ == 1){
        	switch (i){
        		case (0):
        			g.x = -1;
        			g.y = 0;
        			break;
   
        		case (1) :
        			g.x = 1;
        			g.y = 0;
        			break;
        			
        			}
        }
        if ( dim_ == 2 ){
        	switch (i){
        		case (0) :
        			g.x = -1;
        			g.y = -1;
        			break;
        			
        		case (1) :
        			g.x = 1;
        			g.y = 0;
        			break;
        			
        		case (2) :
        			g.x = 0;
        			g.y = 1;
        			break;
        			
        			}
        } 
        
        return g ;
    }

    /****************************************************************/
    /* Implementation of Finite Element functions */
    /****************************************************************/
    void assemble_elementary_matrix(
        const ElementMapping& elt_mapping,
        const ShapeFunctions& reference_functions,
        const Quadrature& quadrature,
        double (*coefficient)(vertex), // pointeur de fonctions
        DenseMatrix& Ke )
    {
        //std::cout << "compute elementary matrix" << '\n';        
        Ke.set_size (reference_functions.nb_functions(),reference_functions.nb_functions());
        for ( int i = 0 ; i < reference_functions.nb_functions()  ; ++i ){        	
        	for ( int j = 0 ; j < reference_functions.nb_functions() ; ++j ){
        		Ke.set(i,j,0.);
        		for ( int q = 0 ; q < quadrature.nb_points() ; ++q ){
        			vertex p_q = quadrature.point(q);
        			double w_q = quadrature.weight(q);
        			DenseMatrix J = elt_mapping.jacobian_matrix(p_q);
        			DenseMatrix J_invert = J.invert_2x2();
        			DenseMatrix J_invert_transpose = J_invert.transpose(); 
        			vec2 gradi = J_invert_transpose.mult_2x2_2(reference_functions.evaluate_grad(i,p_q));
        			vec2 gradj = J_invert_transpose.mult_2x2_2(reference_functions.evaluate_grad(j,p_q));
   			      			       			        			        			
        			Ke.add(i,j,w_q 
        			           * coefficient(elt_mapping.transform(p_q))
        			           * dot (gradi, gradj)
        			           * elt_mapping.jacobian(p_q));        			        
        			}      	          			
        		}	
        	}
        }
    

    void local_to_global_matrix(
        const Mesh& M,
        int t,
        const DenseMatrix& Ke,
        SparseMatrix& K )
    {
        std::cout << "Ke -> K" << '\n';               
        for (int i = 0 ; i < Ke.height() ; ++i){  
        	int i_global = M.get_triangle_vertex_index(t,i);
        	for (int j = 0 ; j < Ke.height() ; ++j){        	      	
        	        int j_global = M.get_triangle_vertex_index(t,j); 
        		K.add(i_global, j_global,Ke.get(i,j));
        	}
        }              
    }
/*
    void assemble_elementary_vector(
        const ElementMapping& elt_mapping,
        const ShapeFunctions& reference_functions,
        const Quadrature& quadrature,
        double (*source)(vertex),
        std::vector< double >& Fe )
    {
        std::cout << "compute elementary vector (source term)" << '\n';
    	Fe.reserve (reference_functions.nb_functions());
        for ( int i = 0 ; i < reference_functions.nb_functions()  ; ++i ){        	
        	for ( int q = 0 ; q < quadrature.nb_points() ; ++q ){
        		vertex p_q = quadrature.point(q);
        		double w_q = quadrature.weight(q);   			      			       			        			        			
        		Fe.push_back(w_q
        		             * reference_functions.evaluate(i,p_q)
        			     * source(elt_mapping.transform(p_q))        			     
        			     * elt_mapping.jacobian(p_q));        			        
        			}      	          			
        		}	
        	}
        }*/

    void assemble_elementary_neumann_vector(
        const ElementMapping& elt_mapping_1D,
        const ShapeFunctions& reference_functions_1D,
        const Quadrature& quadrature_1D,
        double (*neumann)(vertex),
        std::vector< double >& Fe )
    {
        std::cout << "compute elementary vector (neumann condition)" << '\n';
        // TODO
    }

    void local_to_global_vector(
        const Mesh& M,
        bool border,
        int i,
        std::vector< double >& Fe,
        std::vector< double >& F )
    {
        std::cout << "Fe -> F" << '\n';
        // TODO
    }

    void apply_dirichlet_boundary_conditions(
        const Mesh& M,
        const std::vector< bool >& attribute_is_dirichlet, /* size: nb of attributes */
        const std::vector< double >& values, /* size: nb of DOFs */
        SparseMatrix& K,
        std::vector< double >& F )
    {
        std::cout << "apply dirichlet boundary conditions" << '\n';
        // P = 10 000
        // parcours de chauqe bord
        double penalty_coefficient = 10000.;
        for (int i = 0 ; i < M.nb_edges() ; ++i ){
        	if (attribute_is_dirichlet[i]){
        	// parcours chaque noeud 
        		
        	}
        
        }
    }

    void solve_poisson_problem(
            const Mesh& M,
            double (*diffusion_coef)(vertex),
            double (*source_term)(vertex),
            double (*dirichlet_fct)(vertex),
            double (*neumann_fct)(vertex),
            std::vector<double>& solution,
            int verbose )
    {
        std::cout << "solve poisson problem" << '\n';
        // TODO
    }

}
