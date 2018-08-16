#include <string>
#include <vector>
#include <iostream>
#include <CPPALP/LatticePolytope.h>

// h11=491 input string is
// "{{1,0,0,0},{0,1,0,0},{0,0,1,0},{21,28,36,42},{-63,-56,-48,-42}}"

int main(int argc, char* argv[]) {
	std::string verts = argv[1];
	LatticePolytope lp(verts);
	IntMatrix dverts = lp.polar().normal_form();
	
	std::cout << "======= Dual Vertices ========" << std::endl;
	std::cout << intmat_to_string(dverts) << std::endl;

	LatticePolytope dlp = LatticePolytope(dverts);

	for (int i = dlp.dim()-1; i >= 2; i--) {
		std::vector<IntMatrix> fv = dlp.faces(i);
		std::cout << "======= Dimension " << i << " ========" << std::endl;
		for (const IntMatrix& fm : fv) {
			//std::cout << fm << "\n=======" << std::endl;
			// std::cout << intmat_to_string(fm) << std::endl;
			LatticePolytope fp(fm);
			//std::cout << fp.normal_form() << "\n=======" << std::endl;
			std::cout << intmat_to_string(fp.normal_form()) << "\n=======" << std::endl;
			//std::cout << "(" << fp.npoints() << ", " << fp.boundary_points().cols() << ", " << fp.interior_points().cols() << ")" << std::endl;
		}
	}
}