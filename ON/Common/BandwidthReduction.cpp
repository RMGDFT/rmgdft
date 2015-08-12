
#include <boost/config.hpp>
#include <vector>
#include <iostream>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/cuthill_mckee_ordering.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/bandwidth.hpp>

#include "rmgtypedefs.h"
#include "transition.h"

void BandwidthReduction(int num_ions, ION *ions, int *perm_index)
{
    using namespace boost;
    using namespace std;
    typedef adjacency_list<vecS, vecS, undirectedS, 
            property<vertex_color_t, default_color_type,
            property<vertex_degree_t,int> > > Graph;
    typedef graph_traits<Graph>::vertex_descriptor Vertex;
    typedef graph_traits<Graph>::vertices_size_type size_type;

    int ion1, ion2;
    double x, y, z, r, radius = 10.0;
    //Graph G(10);
    Graph G;
    for(ion1 = 0; ion1 < num_ions; ion1++)
    {
        for(ion2 = 0; ion2 < num_ions; ion2++)
        {
            x = fabs(ions[ion1].crds[0] - ions[ion2].crds[0]);
            y = fabs(ions[ion1].crds[1] - ions[ion2].crds[1]);
            z = fabs(ions[ion1].crds[2] - ions[ion2].crds[2]);
            if(x > Rmg_L.a0[0]/2.0) x = Rmg_L.a0[0]-x;
            if(y > Rmg_L.a1[1]/2.0) y = Rmg_L.a1[1]-y;
            if(z > Rmg_L.a2[2]/2.0) z = Rmg_L.a2[2]-z;
            r = sqrt(x*x +y*y + z*z);

            if (r < radius * 2.0)
            {
                add_edge(ion1, ion2, G);
            }
        }
    }

    //  graph_traits<Graph>::vertex_iterator ui, ui_end;

    //  property_map<Graph,vertex_degree_t>::type deg = get(vertex_degree, G);
    //  for (boost::tie(ui, ui_end) = vertices(G); ui != ui_end; ++ui)
    //    deg[*ui] = degree(*ui, G);
    property_map<Graph, vertex_index_t>::type
        index_map = get(vertex_index, G);

    std::cout << index_map[0] << std::endl;
    std::cout << "original bandwidth: " << bandwidth(G) << std::endl;

    std::vector<Vertex> inv_perm(num_vertices(G));
    std::vector<size_type> perm(num_vertices(G));

    {
        //reverse cuthill_mckee_ordering
        Vertex s = vertex(0, G);
        cuthill_mckee_ordering(G, s, inv_perm.rbegin(), get(vertex_color, G),
                make_degree_map(G));

        cout << "Reverse Cuthill-McKee ordering:" << endl;

        for (size_type c = 0; c != inv_perm.size(); ++c)
            perm[index_map[inv_perm[c]]] = c;
        std::cout << "  bandwidth: " 
            << bandwidth(G, make_iterator_property_map(&perm[0], index_map, perm[0]))
            << std::endl;
    }

    unsigned int i;
    for (i = 0; i < inv_perm.size(); i++) perm_index[i] = inv_perm[i];

    #if 0
    for(i = 0; i < num_ions; i++)
    {
        ion1 = inv_perm[i];
    for(int j = 0; j < num_ions; j++)
        {
        ion2 = inv_perm[j];
            x = fabs(ions[ion1].crds[0] - ions[ion2].crds[0]);
            y = fabs(ions[ion1].crds[1] - ions[ion2].crds[1]);
            z = fabs(ions[ion1].crds[2] - ions[ion2].crds[2]);
            //if(x > Rmg_L.a0[0]/2.0) x = Rmg_L.a0[0]-x;
            if(y > Rmg_L.a1[1]/2.0) y = Rmg_L.a1[1]-y;
            if(z > Rmg_L.a2[2]/2.0) z = Rmg_L.a2[2]-z;
            r = sqrt(x*x +y*y + z*z);

            if (r < radius * 2.0)
                printf("\n  %d  %d  abababab", i, j);
        }
    }
    #endif
}
