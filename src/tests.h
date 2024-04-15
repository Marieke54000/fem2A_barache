#pragma once

#include "mesh.h"
#include "fem.h"
#include "solver.h"

#include <assert.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <stdlib.h>

namespace FEM2A {
    namespace Tests {

        bool test_load_mesh()
        {
            Mesh mesh;
            mesh.load("data/square.mesh");

            std::cout << "Vertices <x> <y> <att>" << std::endl;
            for( int v = 0; v < mesh.nb_vertices(); v++ ) {
                std::cout << mesh.get_vertex(v).x << " " << mesh.get_vertex(v).y
                    << " " << mesh.get_vertex_attribute(v) << std::endl;
            }

            std::cout << "Edges <v0> <v1> <att>" << std::endl ;
            for( int ed = 0; ed < mesh.nb_edges(); ed++ ) {
                std::cout << mesh.get_edge_vertex_index(ed, 0) << " "
                    << mesh.get_edge_vertex_index(ed, 1) << " "
                    << mesh.get_edge_attribute(ed) << std::endl;
            }

            std::cout << "Triangles <v0> <v1> <v2> <att>" << std::endl ;
            for( int tr = 0; tr < mesh.nb_triangles(); tr++ ) {
                std::cout << mesh.get_triangle_vertex_index(tr, 0) << " "
                    << mesh.get_triangle_vertex_index(tr, 1) << " "
                    << mesh.get_triangle_vertex_index(tr, 2) << " "
                    << mesh.get_triangle_attribute(tr) << std::endl;
            }

            return true;
        }

        bool test_load_save_mesh()
        {
            Mesh mesh;
            mesh.load("data/geothermie_4.mesh");
            mesh.save("data/geothermie_4.mesh");
            return true;
        }
        
        /**
        Ce test verifie que la somme des poids pour toute fonction f=1 vaut 1/2
        **/
        
        bool test_quadrature( int order)  
        {
            Quadrature quad = Quadrature::get_quadrature(order);
            std :: cout << quad.nb_points() << std :: endl;
            double sum =0;
            for (int i =0; i < quad.nb_points(); ++i ){
            	std::cout << quad.point(i).x << " " // point a été défini commme une structure vertex (dans mesh.h)
            		<< quad.point(i).y << std::endl;
            	std::cout << quad.weight(i) << std::endl;
            	sum = sum + quad.weight(i);
            	}
            std::cout << sum << std::endl;
            return true;
        }
       
       bool test_element_mapping( bool border, int i )
       {
       	    Mesh mesh;
       	    mesh.load("data/square.mesh");
       	    ElementMapping map = ElementMapping ( mesh,  border, i );
       	    vertex x_r ; 
       	    
       	    vertex r = map.transform(x_r);
       	    std :: cout << r.x << " " << r.y << std :: endl ;       	    
       	    DenseMatrix J_r = map.jacobian_matrix (x_r);
       	    std :: cout << J_r.get (0,0) << " " << J_r.get (0,1) << '\n';
       	    std :: cout << J_r.get (1,0) << " " << J_r.get (1,1) << '\n';
       	    std :: cout << map.jacobian (x_r) << '\n';
       	    return true ;
       		
    	}
    	
    	bool test_shape_functions (int dim, int order)
    	{
            ShapeFunctions fct = ShapeFunctions (dim , order);
            std :: cout << "il y a " << fct.nb_functions() << " fonctions " << std :: endl ;           
            int i;
            vertex x_r ;
            x_r.x = 0.2 ;
       	    x_r.y = 0.4 ;
            std :: cout << " la valeur de la fonction de forme en ce point est : " << fct.evaluate(0, x_r) << std ::endl;
            return true;
        }
   
    }
}
