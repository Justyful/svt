#include "inmost.h"
#include <stdio.h>

double city_x = 109.42622950819671;
double city_y = -6.097560975609753;

using namespace INMOST;
using namespace std;

void make_vec_to_center_tag(Mesh *m)
{
	Tag vec_city = m->CreateTag("Vector from cell center to Serpukhov", DATA_REAL, CELL, NONE, 3);
	for (Mesh::iteratorCell icell = m->BeginCell(); icell != m->EndCell(); icell++) {
		Cell cell = icell->getAsCell();
		double *coords = new double[3];
		cell.Barycenter(coords);
		cell.RealArray(vec_city)[0] = city_x - coords[0];
		cell.RealArray(vec_city)[1] = city_y - coords[1];
		cell.RealArray(vec_city)[2] = 0;
	}
}


double dist_between_nodes(Node node1, Node node2) {
	double a = node1.Coords()[0] - node2.Coords()[0];
	double b = node1.Coords()[1] - node2.Coords()[1];
	return sqrt(a * a + b * b);
}


double mesh_diam(Mesh* m)
{
	double max_diam = 0;
	for (Mesh::iteratorCell icell = m->BeginCell(); icell != m->EndCell(); icell++) {
		Cell cell = icell->getAsCell();
		ElementArray<Node> nodes = cell.getNodes();

		double dist01 = dist_between_nodes(nodes[0], nodes[1]);
		double dist02 = dist_between_nodes(nodes[0], nodes[2]);
		double dist12 = dist_between_nodes(nodes[1], nodes[2]);

		double cell_diam = max(max(dist01, dist02), dist12);
		if (cell_diam > max_diam) {
			max_diam = cell_diam;
		}
	}
	return max_diam;
}

void make_cells_count_tag(Mesh* m)
{
	Tag cells_cnt = m->CreateTag("Number of adjacent cells", DATA_INTEGER, NODE, NONE, 1);
	for (Mesh::iteratorNode inode = m->BeginNode(); inode != m->EndNode(); inode++) {
		Node node = inode->getAsNode();
		node.Integer(cells_cnt) = node.getCells().size();
	}
}

int main(int argc, char *argv[])
{
	if (argc != 2) {
		cout << "Usage: mesh <mesh.vtk>" << endl;
		return 1;
	}
    Mesh *m = new Mesh;
	m->Load(argv[1]);
	make_cells_count_tag(m);
	m->Save("res.vtk");
	delete m;
	cout << "Success! \n";
	return 0;
}