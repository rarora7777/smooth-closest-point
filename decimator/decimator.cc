#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <sstream>
#include <string>
#include <memory>
#include <map>


#define OM_STATIC_BUILD
#define _USE_MATH_DEFINES
#include <OpenMesh/Core/IO/MeshIO.hh>
//#include <OpenMesh/Core/IO/writer/OFFWriter.hh>
//#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Tools/Decimater/DecimaterT.hh>
#include <OpenMesh/Tools/Decimater/ModAspectRatioT.hh>
#include <OpenMesh/Tools/Decimater/ModEdgeLengthT.hh>
#include <OpenMesh/Tools/Decimater/ModHausdorffT.hh>
#include <OpenMesh/Tools/Decimater/ModNormalDeviationT.hh>
#include <OpenMesh/Tools/Decimater/ModNormalFlippingT.hh>
#include <OpenMesh/Tools/Decimater/ModQuadricT.hh>
#include <OpenMesh/Tools/Decimater/ModProgMeshT.hh>
#include <OpenMesh/Tools/Decimater/ModIndependentSetsT.hh>
#include <OpenMesh/Tools/Decimater/ModRoundnessT.hh>

//typedef OpenMesh::DefaultTraits MyTraits;
//typedef OpenMesh::PolyMesh_ArrayKernelT<>  MyMesh;
typedef OpenMesh::TriMesh_ArrayKernelT<> MyMesh;
typedef OpenMesh::Decimater::DecimaterT<MyMesh> Decimater;
typedef OpenMesh::Decimater::ModQuadricT<MyMesh>::Handle HModQuadric;
typedef OpenMesh::Decimater::ModNormalDeviationT<MyMesh>::Handle HModNormalDeviation;
typedef OpenMesh::Decimater::ModEdgeLengthT<MyMesh>::Handle HModEdgeLengthT;

#include <OpenMesh/Core/Utils/vector_cast.hh>

int decimate (char* meshfile, unsigned stop_at_vertices, bool normaldeviation, bool edgelength)
{
    MyMesh mesh;

    if (!OpenMesh::IO::read_mesh (mesh, meshfile))
    {
	std::cerr << "Error loading mesh from file " << meshfile << std::endl;
	exit (-1);
    }

    if ( !mesh.has_face_normals() )
      mesh.request_face_normals();

    mesh.update_face_normals();

    if ( !mesh.has_vertex_normals() )
      mesh.request_vertex_normals();

    mesh.update_vertex_normals();

	OpenMesh::VPropHandleT<int> ids;
  	mesh.add_property(ids);
	{  	
		unsigned i = 0;
		MyMesh::VertexIter v_it, v_end(mesh.vertices_end());
    	for (v_it=mesh.vertices_begin(); v_it!=v_end; ++v_it)    
    	{
    		mesh.property(ids,*v_it) = i++;
    	}
 	} 	

    Decimater decimater (mesh);
    
    HModQuadric hModQuadric;
    HModNormalDeviation hModNormalDeviation;
    HModEdgeLengthT hModEdgeLengthT;

    
    if (normaldeviation)
    {
    	decimater.add(hModNormalDeviation);
    	decimater.module( hModNormalDeviation ).set_binary( false );
    }	
    else if (edgelength)
    {
    	decimater.add(hModEdgeLengthT);
    	decimater.module( hModEdgeLengthT ).set_binary( false );
    }
    else
    {
    	decimater.add (hModQuadric);
    }
    	
    bool rc = decimater.initialize ();
    
    if (!rc)
    {
         std::cerr << "  initializing failed!" << std::endl;
         decimater.info( std::cerr );
         return false;
    }

 
    int level = 0;
    unsigned prev_num_vertices = 0;

	int ncoll = decimater.decimate_to ((int)(stop_at_vertices));

	mesh.garbage_collection();
	
    std::stringstream ss2;
    ss2 << meshfile << "-decimated.map";

    FILE *fp = fopen (ss2.str().c_str(), "w");

	{
		MyMesh::VertexIter v_it, v_end(mesh.vertices_end());
    	for (v_it=mesh.vertices_begin(); v_it!=v_end; ++v_it)    
    	{
   			fprintf (fp, "%d ", mesh.property(ids,*v_it));
    	}
    	fprintf (fp, "\n");
    	fclose (fp);
    }



    std::stringstream ss;
    ss << meshfile << "-decimated.off";
    if (!OpenMesh::IO::write_mesh (mesh, ss.str().c_str()))
    {
		std::cerr << "Cannot write mesh to file '" << ss.str().c_str() << "'" << std::endl;
		exit (-1);
    }


    return level;
}

void print_help ()
{
    fprintf (stderr, "Usage: decimator -m input_file [-d decimation_factor] [-s stop_at_vertices] [-n] [-l]\n"
             "\t\t -m input_file: Original mesh in off format.\n"
             "\t\t -s stop_at_vertices: cut until this number of vertices is reached, then stop. Integer value greater than 10, default 100.\n"
             "\t\t -n Optimize for Normal Deviation .\n"
             "\t\t -l Optimize for Edge Length.\n"
	     "\t\t -h: print this help and exit.\n");
}

char* meshfile = 0;
int stop_at_vertices = 100;
bool normaldeviation = false;
bool edgelength = false;

void app_init(int argc, char *argv[])
{
    for (argc--, argv++; argc--; argv++)
    {
        if( (*argv)[0] == '-')
        {
            switch( (*argv)[1] )
            {
            case 'm':
                meshfile = (*++argv);
                argc--;
                break;

            case 's':
                stop_at_vertices = atoi (*++argv);
                argc--;
                break;

            case 'n':
                normaldeviation = true;
                break;

            case 'l':
                edgelength = true;
                break;

            case 'h':
                print_help();
                exit(0);

            default:
                goto die;
            }
        }
    }

    if (!meshfile)
	goto die;

    char *p;
    
    if (stop_at_vertices <= 10)
	goto die;
    
    return;

 die:
    print_help();
    exit(1);
}

int main(int argc, char *argv[])
{
    int levels = 0;
    app_init (argc, argv);

    levels = decimate (meshfile, stop_at_vertices,normaldeviation,edgelength);

    exit (0);
}
