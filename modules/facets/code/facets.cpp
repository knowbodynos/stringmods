/*

 ____ ______ ___________  |  | ______  
_/ ___\\____ \\____ \__  \ |  | \____ \ 
\  \___|  |_> >  |_> > __ \|  |_|  |_> >
 \___  >   __/|   __(____  /____/   __/ 
     \/|__|   |__|       \/     |__|   


    ---------------------------------
	|	 Created by Jon Carifio		|
	| 		January 2018			|
    -------------------------------*/

/*
 The purpose of this file is to provide a C++ wrapper for several of the primary PALP 
functions involving lattice polyhedra. The currently implemented features are

-Polytope dimension
-Number of vertices
-Number of integral points
-Number of facets
-Check reflexivity
-Polar dual (if reflexive)
-Hodge numbers (if reflexive)
-Normal form (sometimes this calculation is weirdly slow)
-f-vector
-Faces of each dimension

as well as some helper functions (i.e. for using with the database)


The PALP functions that are needed for the current version can all be found in the "Global6.h"
header file. In order to distinguish the native PALP functions from the ones create in this file, 
they are wrapped in the namespace PALP.
*/

#include <exception>
#include <vector>
#include <string>
#include <iostream>
#include <bitset>
#include <algorithm>
#include <stdio.h>
#include "Eigen/Core"
#include "Eigen/LU"
#include "rapidjson/document.h"

using namespace rapidjson;
using namespace std;

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> PolyMatrix;

namespace PALP {
	#include "GlobalP.h"
}

// Forward declarations

// Printing functions (for testing purposes)
void print_ppl(PALP::PolyPointList* ppl, int amb);
void print_vnl(PALP::VertexNumList* vnl);
void print_eqlist(PALP::EqList* eqlist);
PolyMatrix vector_to_polymat(vector<vector<int>> v);


// Exception for methods that require reflexivity
class NotReflexiveException : public exception {
	virtual const char* what() const throw()
  {
    return "This function is implemented only for reflexive polytopes!";
	}
};

// Exception if the user inputs a point with the wrong ambient dimension
class WrongDimensionException : public exception {
	virtual const char* what() const throw()
  {
    return "The dimension of the input point does not match that of the polytope's ambient space.";
	}
};

// Exception if the polytope needs to be full dimension
class FullDimensionException : public exception {
	virtual const char* what() const throw()
  {
    return "Only implemented for full-dimensional polytopes.";
	}
};


class LatticePolytope {

	public:
		// Constructors
		LatticePolytope(PolyMatrix pmat);
		LatticePolytope(vector<vector<int>>v ) : LatticePolytope(vector_to_polymat(v)) {}

		// Destructor
		// Currently there isn't anything special that needs to be done in the destructor
		// Any 'extra' storage created in the constructor is destroyed there
		~LatticePolytope() {}

		// Obtain member values
		int ambient_dim() {return amb;}
		int npoints() {return np;}
		int nvertices() {return nv;}
		int nfacets() {return ne;}
		int dim() {return dmn;}
		bool is_reflexive() {return isref;}
		

		// Obtain properties
		PolyMatrix vertices();
		PolyMatrix integral_points();
		PolyMatrix vertex_edge_pairing_matrix();
		vector<PolyMatrix> faces(const int& fdim);
		vector<int> f_vector();
		vector<int> hodge();


		// Obtain PALP member values
		PALP::PolyPointList get_plist() {return plist;}
		PALP::VertexNumList get_vnlist() {return vnlist;}
		PALP::EqList get_eqlist() {return eqlist;}
		//PALP::PairMat get_vepm() {return &vepm;}

		// Print PALP member values
		void print_equations();

		// Calculations
		LatticePolytope polar();
		PolyMatrix normal_form();
		PolyMatrix normal_form2(); //don't use
		bool boundary_contains(const vector<int>& pt);
		bool interior_contains(const vector<int>& pt);
		bool p_boundary_contains(const vector<int>& pt, const int& p);
		bool p_interior_contains(const vector<int>& pt, const int& p);
		PolyMatrix interior_points();
		PolyMatrix boundary_points();
		PolyMatrix p_interior_points(const int& p);
		PolyMatrix p_boundary_points(const int& p);

		// Get the vertices in database format
		string vertex_to_string(const int& n);
		string vertex_string();
		

	private:

		// Max number for bitstring
		const static int BIT_MAX = 64;
		
		// Data members
		int amb;
		int dmn;
		int nv;
		int np;
		int ne;
		int isref;

		// PALP data members (used for calling PALP functions)
		PALP::PolyPointList plist;
		PALP::VertexNumList vnlist;
		PALP::EqList eqlist;
		PALP::PairMat vepm;
		PALP::FaceInfo finfo;

		// Internal calculations
		vector<int> INCI_to_indices(const PALP::INCI& inc);

};


// Constructor
LatticePolytope::LatticePolytope(PolyMatrix pmat)

	//Trivial data members
	:amb(pmat.rows()), nv(pmat.cols()),

	vepm{{0}}
	{

		// Get the dimension
		PolyMatrix dmat(pmat.rows(), pmat.cols()-1);
		for (int i = 1; i < pmat.cols(); i++) {
			for (int j = 0; j < pmat.rows(); j++) {
				dmat(j,i-1) = pmat(j,i) - pmat(j,0);
			}
		}

		Eigen::FullPivLU<PolyMatrix> lu = dmat.fullPivLu();
		dmn = lu.rank();

		// Construct the PolyPointList pointer
		PALP::PolyPointList* vlist = new PALP::PolyPointList;
		vlist->n = dmn;
		vlist->np = pmat.cols();
		for (int i = 0; i < vlist->np; i++) {
			for (int j = 0; j < amb; j++) {
				vlist->x[i][j] = (long) pmat(j,i);
				// cout << pplist->x[i][j] << "\t";
			}
			// cout << endl;
		}
		//cout << vlist->n << ", " << vlist->np << endl;
		
		// Construct the VertexNumList pointer
		// Initially, we set this to include every point - the correct indices are found later
		/*PALP::VertexNumList* vnl = new PALP::VertexNumList;
		vnl->nv = nv;
		for (int i=0; i < nv; i++) {
			vnl->v[i] = i;
			// cout << vnl->v[i] << endl;
		}
		vnlist = *vnl;*/
		
		/*
		PolyMatrix q(pmat.rows(), pmat.cols());
		PolyMatrix r(pmat.rows(), 1);
		for (int i = 0; i < pmat.cols(); i++) {
			r = (dmat.transpose().fullPivLu().matrixL()*dmat).transpose()*(dmat.transpose().fullPivLu().matrixL()*dmat).transpose().fullPivLu().solve(pmat.col(i)-pmat.col(0));
			for (int j = 0; j < pmat.rows(); j++) {
				q(j,i) = r(j,0)+pmat(j,0);
			}
		}
		cout << (dmat.transpose().fullPivLu().matrixL()*dmat).transpose() << endl << endl;
		cout << q << endl;
		*/

		// Construct the Equation pointer
		// Determine reflexivity on the way
		PALP::EqList* eql = new PALP::EqList;
		PALP::VertexNumList* vnl = new PALP::VertexNumList;
		isref = PALP::Ref_Check(vlist, vnl, eql);

		// If the polytope isn't reflexive, we need to calculate the equations
		if (!isref) {
			PALP::Find_Equations(vlist, vnl, eql);
		}

		eqlist = *eql;
		vnlist = *vnl;

		// Construct the PALP vertex pairing matrix
		PALP::Make_VEPM(vlist, vnl, &eqlist, vepm);

		// Construct the full list of points
		//cout << "About to complete polytope!" << endl;
		PALP::PolyPointList* pplist = new PALP::PolyPointList;
		pplist->n = dmn;
		pplist->np = nv;
		for (int i = 0; i < nv; i++) {
			for (int j = 0; j < amb; j++) {
				pplist->x[i][j] = vlist->x[i][j];
			}
		}
		PALP::Complete_Poly(vepm, &eqlist, nv, pplist);
		plist = *pplist;
		
		// Set the number of points and equations
		np = plist.np;
		ne = eqlist.ne;

		// Trim vlist, if necessary
		vector<int> vcheck;
		PALP::VertexNumList* vnl2 = new PALP::VertexNumList;
		for (int i = 0; i < vnlist.nv; i++) {
			int val = vnlist.v[i];
			if (!(find(vcheck.begin(), vcheck.end(), val) != vcheck.end())) {
				vcheck.push_back(val);
			}
		}
		// sort(vcheck.begin(),vcheck.end());

		vnl2->nv = vcheck.size();
		for (int i = 0; i < vcheck.size(); i++) {
			vnl2->v[i] = vcheck[i];
		}
		vnlist = *vnl2;

		// Final nv assignment
		nv = vnlist.nv;

		// Get the incidence face info
		//cout << "Time for the face info!" << endl;
		PALP::FaceInfo* finf = new PALP::FaceInfo;
		PALP::Make_Incidence(&plist, &vnlist, &eqlist, finf);
		finfo = *finf;

		// Clear any necessary memory
		delete vnl;
		delete eql;
		delete vlist;
		delete pplist;
		delete vnl2;
		delete finf;

	}


// Class member functions

LatticePolytope LatticePolytope::polar() {

	// Check for reflexivity
	if (!is_reflexive()) {throw NotReflexiveException();}

	PALP::PolyPointList* dpl = new PALP::PolyPointList;
	PALP::Make_Dual_Poly(&plist, &vnlist, &eqlist, dpl);

	PolyMatrix dpm(amb,dpl->np);
	for (int i = 0; i < dpl->np; i++) {
		for (int j = 0; j < amb; j++) {
			dpm(j,i) = (double) dpl->x[i][j];
		}
	}
	delete dpl;
	return LatticePolytope(dpm);
}

PolyMatrix LatticePolytope::normal_form2() {

	// Check for reflexivity
	//if (!is_reflexive()) {throw NotReflexiveException();}

	PolyMatrix nf(4, nv);
	long pNF[POLY_Dmax][VERT_Nmax];
	int sn = 0;
	int* symnum = &sn;
	int V_perm[SYM_Nmax][VERT_Nmax];
	int t = 0;
	int S = 0;
	int N = 0;
	//PALP::Make_Poly_NF(plist, vnlist, eqlist, pNF);
	PALP::Make_Poly_Sym_NF(&plist, &vnlist, &eqlist, symnum, V_perm, pNF, t, S, N);
	for (int i = 0; i < nv;i ++) {
		for (int j = 0; j < amb; j++) {
			nf(j,i) = pNF[j][i];
		}
	}
	return nf;
}

PolyMatrix LatticePolytope::normal_form() {

	// Check if the polytope is full-dimensional
	if (amb != dmn) {throw FullDimensionException();}

	PolyMatrix nf(4, nv);
	long pNF[POLY_Dmax][VERT_Nmax];
	/*PALP::EqList equations = *eqlist;
	PALP::PolyPointList ppl = *plist;
	PALP::VertexNumList vnl = *vnlist;*/
	PALP::Make_Poly_NF(&plist, &vnlist, &eqlist, pNF);
	/*plist = &ppl;
	vnlist = &vnl;
	eqlist = &equations;*/
	for (int i = 0; i < nv; i++) {
		for (int j = 0; j < amb; j++) {
			nf(j,i) = pNF[j][i];
		}
	}
	return nf;
	
}

PolyMatrix LatticePolytope::vertices() {
	PolyMatrix vm(amb, nv);
	for (int i = 0; i < nv; i++) {
		int index = vnlist.v[i];
		for (int j = 0; j < amb; j++) {
			vm(j,i) = plist.x[index][j];
		}
	}
	return vm;
}

PolyMatrix LatticePolytope::integral_points() {
	PolyMatrix ipm(amb, np);
	for (int i = 0; i < np; i++) {
		for (int j = 0; j < amb; j++) {
			ipm(j,i) = plist.x[i][j];
		}
	}
	return ipm;
}

PolyMatrix LatticePolytope::vertex_edge_pairing_matrix() {
	PolyMatrix vpm(ne, nv);
	for (int i = 0; i < vpm.rows(); i++) {
		for (int j = 0; j < vpm.cols(); j ++) {
			vpm(i,j) = vepm[i][j];
		}
	}
	return vpm;
}

vector<int> LatticePolytope::INCI_to_indices(const PALP::INCI& inc) {
	// This function converts an INCI object (an incidence bit pattern) to a subset of vnlist
	string s = bitset<BIT_MAX>((unsigned long long)inc).to_string();
	vector<int> indices;
	for (int i=0; i < s.length(); i++) {
		if (s.substr(i,1) == "1") {
			//indices.push_back(BIT_MAX-i-1);
			indices.push_back(nv-BIT_MAX+i);
		}
	}
	sort(indices.begin(), indices.end());
	return indices;
}

vector<PolyMatrix> LatticePolytope::faces(const int& fdim) {
	vector<PolyMatrix> fvec;
	int nfaces = finfo.nf[fdim];
	for (int i = 0; i < nfaces; i++) {
		PALP::INCI bitp = finfo.v[fdim][i];
		vector<int> vec_inds = INCI_to_indices(bitp);
		PolyMatrix fm(amb, vec_inds.size());
		for (int j = 0; j < vec_inds.size(); j++) {
			for (int k = 0; k < amb; k++) {
				fm(k,j) = plist.x[vnlist.v[vec_inds[j]]][k];
			}
		}
		fvec.push_back(fm);
	}
	return fvec;
}

vector<int> LatticePolytope::f_vector() {
	vector<int> fvec;
	fvec.push_back(1);
	for (int i=0; i < dmn; i++) {
		fvec.push_back(finfo.nf[i]);
	}
	fvec.push_back(1);
	return fvec;
}

vector<int> LatticePolytope::hodge() {

	// Check for reflexivity
	if (!is_reflexive()) {throw NotReflexiveException();}

	// Create a container for the Hodge numbers
	PALP::BaHo* baho = new PALP::BaHo;

	// Get the dual polytope info
	PALP::PolyPointList* dpl = new PALP::PolyPointList;
	PALP::Make_Dual_Poly(&plist, &vnlist, &eqlist, dpl);
	PALP::RC_Calc_BaHo(&plist, &vnlist, &eqlist, dpl, baho);

	// Create the return vector
	vector<int> hnums;
	hnums.push_back(baho->h1[1]);
	hnums.push_back(baho->h1[2]);

	// Free space and return
	//free(baho->h1);
	delete dpl;
	delete baho;
	return hnums;

}

bool LatticePolytope::boundary_contains(const vector<int>& pt) {
	if (pt.size() != amb) {throw WrongDimensionException();}

	/* 
	For each equation that defines the edges of the polytope, we see if the point satisfies
	at least one. If it does, it lies on the boundary of the polytope. Otherwise it doesn't.
	*/
	for (int j = 0; j < eqlist.ne; j++) {
		int sp = 0;
		PALP::Equation eq = eqlist.e[j];
		for (int i = 0; i < amb; i++) {
			sp += pt[i]*eq.a[i];
		}

		if (sp == -eq.c) {
			return true;
		}
	}
	return false;
}


bool LatticePolytope::interior_contains(const vector<int>& pt) {
	if (pt.size() != amb) {throw WrongDimensionException();}
	if (amb != dmn) {throw FullDimensionException();}

	/* 
	For each hyperplane that defines the edges of the polytope, we see if the point lies inside 
	of it. If it is inside all of them, it lies within the interior. Otherwise it doesn't.
	Because of this, degenerate (not full-dimensional) cases aren't yet supported.
	*/
	for (int j = 0; j < eqlist.ne; j++) {
		int sp = 0;
		PALP::Equation eq = eqlist.e[j];
		for (int i = 0; i < amb; i++) {
			sp += pt[i]*eq.a[i];
		}

		if (sp <= -eq.c) {
			return false;
		}
	}
	return true;
}


bool LatticePolytope::p_boundary_contains(const vector<int>& pt, const int& p) {
	if (pt.size() != amb) {throw WrongDimensionException();}

	/* 
	For each equation that defines the edges of the polytope, we see if the point satisfies
	at least one. If it does, it lies on the boundary of the polytope. Otherwise it doesn't.
	*/
	int eq_count = 0;
	for (int j = 0; j < eqlist.ne; j++) {
		PALP::Equation eq = eqlist.e[j];
		int sp = 0;
		for (int i = 0; i < amb; i++) {
			sp += pt[i]*eq.a[i];
		}

		if (sp == -eq.c) {
			eq_count++;
			if (eq_count == amb - p + 1) {
				return true;
			}
		}
	}
	return false;
}


bool LatticePolytope::p_interior_contains(const vector<int>& pt, const int& p) {
	if (pt.size() != amb) {throw WrongDimensionException();}
	if (amb != dmn) {throw FullDimensionException();}

	/* 
	For each hyperplane that defines the edges of the polytope, we see if the point lies inside 
	of it. If it is inside all of them, it lies within the interior. Otherwise it doesn't.
	Because of this, degenerate (not full-dimensional) cases aren't yet supported.
	*/
	int eq_count = 0;
	for (int j = 0; j < eqlist.ne; j++) {
		PALP::Equation eq = eqlist.e[j];
		int sp = 0;
		for (int i = 0; i < amb; i++) {
			sp += pt[i]*eq.a[i];
		}

		if (sp == -eq.c) {
			eq_count++;
		}
	}
	return eq_count == amb - p;
}


PolyMatrix LatticePolytope::boundary_points() {

	// Determine which points lie on the boundary
	vector<int> bd_indices;
	for (int i = 0; i < np; i++) {
		vector<int> pt;
		for (int j = 0; j < amb; j++) {
			pt.push_back(plist.x[i][j]);
		}
		if (boundary_contains(pt)) {bd_indices.push_back(i);}
	}

	int nbp = bd_indices.size(); // The number of boundary points
	PolyMatrix bdp(amb,nbp);
	for (int j = 0; j < nbp; j++) {
		for (int k = 0; k < amb; k++) {
			bdp(k,j) = plist.x[bd_indices[j]][k];
		}
	}
	return bdp;
}


PolyMatrix LatticePolytope::interior_points() {

	// Determine which points lie in the interior
	vector<int> intr_indices;
	for (int i = 0; i < np; i++) {
		vector<int> pt;
		for (int j = 0; j < amb; j++) {
			pt.push_back(plist.x[i][j]);
		}
		if (!boundary_contains(pt)) {intr_indices.push_back(i);}
	}

	int nintrp = intr_indices.size(); // The number of interior points
	PolyMatrix intrp(amb,nintrp);
	for (int j = 0; j < nintrp; j++) {
		for (int k = 0; k < amb; k++) {
			intrp(k,j) = plist.x[intr_indices[j]][k];
		}
	}
	return intrp;
}


PolyMatrix LatticePolytope::p_boundary_points(const int& p) {

	// Determine which points lie on the boundary
	vector<int> pbd_indices;
	for (int i = 0; i < np; i++) {
		vector<int> pt;
		for (int j = 0; j < amb; j++) {
			pt.push_back(plist.x[i][j]);
		}
		if (p_boundary_contains(pt, p)) {pbd_indices.push_back(i);}
	}

	int npbp = pbd_indices.size(); // The number of boundary points
	PolyMatrix pbdp(amb,npbp);
	for (int j = 0; j < npbp; j++) {
		for (int k = 0; k < amb; k++) {
			pbdp(k,j) = plist.x[pbd_indices[j]][k];
		}
	}
	return pbdp;
}


PolyMatrix LatticePolytope::p_interior_points(const int& p) {

	// Determine which points lie in the interior
	vector<int> pintr_indices;
	for (int i = 0; i < np; i++) {
		vector<int> pt;
		for (int j = 0; j < amb; j++) {
			pt.push_back(plist.x[i][j]);
		}
		if (p_interior_contains(pt, p)) {pintr_indices.push_back(i);}
	}

	int npintrp = pintr_indices.size(); // The number of interior points
	PolyMatrix pintrp(amb,npintrp);
	for (int j = 0; j < npintrp; j++) {
		for (int k = 0; k < amb; k++) {
			pintrp(k,j) = plist.x[pintr_indices[j]][k];
		}
	}
	return pintrp;
}


string LatticePolytope::vertex_to_string(const int& n) {
	vector<int> v;

	// Get the coordinates of the nth vertex
	int vertex_ind = vnlist.v[n];
	string vstr = "{";
	for (int i = 0; i < amb-1; i++) {
		vstr += to_string(plist.x[vertex_ind][i]);
		vstr += ",";
	}
	vstr += to_string(plist.x[vertex_ind][amb-1]);
	vstr += "}";
	return vstr;
}

string LatticePolytope::vertex_string() {
	string vstr = "{";
	for (int i = 0; i < nv-1; i++) {
		vstr += vertex_to_string(i);
		vstr += ",";
	}
	vstr += vertex_to_string(nv-1);
	vstr += "}";
	return vstr;
}

void LatticePolytope::print_equations() {
	cout << eqlist.ne << endl;
	for (int i = 0; i < eqlist.ne; i++) {
		PALP::Equation eq = eqlist.e[i];
		cout << eq.c << endl;
		for (long w : eq.a) {
			cout << w << "\t";
		}
		cout << endl;
	}

}




// Other (helper) functions

string col_to_string(const PolyMatrix&pm, const int& i) {
	vector<int> v;
	string vstr = "{";
	for (int j = 0; j < pm.rows()-1; j++) {
		vstr += to_string((int)pm(j,i));
		vstr += ",";
	}
	vstr += to_string((int)pm(pm.rows()-1,i));
	vstr += "}";
	return vstr;
}

string polymat_to_string(const PolyMatrix& pm) {
	string pmstr = "{";
	for (int i = 0; i < pm.cols()-1; i++) {
		pmstr += col_to_string(pm, i);
		pmstr += ",";
	}
	pmstr += col_to_string(pm, pm.cols()-1);
	pmstr += "}";
	return pmstr;
}

vector<string> split(string s, const string& dlm) {
	// This function splits the given string using the given delimiter
	// Returns a vector containing the pieces of the split string

	// Container to store split pieces
	vector<string> pieces;

	// Length of the delimiter
	int dlm_len = dlm.length();

	// While the delimiter still exists in the string, keep splitting off pieces
	size_t pos = 0;
	string piece;
	while ((pos=s.find(dlm)) != string::npos) {
		piece = s.substr(0, pos);
		pieces.push_back(piece);
		s.erase(0, pos + dlm_len);
	}

	pieces.push_back(s);
	return pieces;

}


vector<int> string_to_int_vector(const vector<string>& svec) {
	// Converts a vector of integers in string for to a vector of ints
	vector<int> intvec;
	for (int i = 0; i < svec.size(); ++i) {
		string s = svec[i];
		int n = stoi(s);
		intvec.push_back(n);
	}
	return intvec;
}

vector<vector<int>> string_to_vertex_list(string s, const string& sep="},{") {
	s.erase(0,2);
	s.erase(s.length()-2,2);
	vector<string> vertices = split(s, sep);
	vector<vector<int>> vertex_list;
	for (string vs : vertices) {
		vector<string> pieces = split(vs, ",");
		vector<int> coords = string_to_int_vector(pieces);
		vertex_list.push_back(coords);
	}
	return vertex_list;
}

PolyMatrix vector_to_polymat(const vector<vector<int>>& v) {
	int npts = v.size();
	int dim = v[0].size();
	PolyMatrix pm(dim, npts);
	for (int i = 0; i < npts; i++) {
		for (int j = 0; j < dim; j++) {
			pm(j,i) = (double) v[i][j];
		}
	}
	return pm;
}

vector<vector<int>> polymat_to_vector(const PolyMatrix& pm) {
	int npts = pm.cols();
	int dim = pm.rows();
	vector<vector<int>> v;
	for (int i = 0; i < npts; i++) {
		vector<int> vv;
		for (int j = 0; j < dim; j++) {
			vv.push_back((int)pm(j,i));
		v.push_back(vv);
		}
	}
	return v;
}

PolyMatrix add_col(const PolyMatrix& pm, const PolyMatrix& col) {
	PolyMatrix new_pm = pm;
	new_pm.conservativeResize(new_pm.rows(), new_pm.cols()+1);
	new_pm.col(new_pm.cols()-1) = col;
	return new_pm;
}

PolyMatrix remove_col(const PolyMatrix& pm, const PolyMatrix& col) {
	vector<int> col_inds;
	for (int i=0; i < pm.cols(); i++)
	{
	    if (pm.col(i) != col) col_inds.push_back(i);
	}
	PolyMatrix temp(pm.rows(), col_inds.size());
	for (int i=0; i < temp.cols(); i++)
	{
		for (int j=0; j < pm.rows(); j++) {
		    temp(j,i) = pm(j,col_inds[i]);
		}
	}
	return temp;
}


// Functions for printing PALP classes

void print_ppl(PALP::PolyPointList* ppl, int amb) {
	for (int i = 0; i < ppl->np; i++) {
		for (int j = 0; j < amb; j++) {
			cout << ppl->x[i][j] << "\t";
		}
		cout << endl;
	}
}

void print_vnl(PALP::VertexNumList* vnl) {
	for (int i = 0; i < vnl->nv; i++) {
		cout << vnl->v[i] << "\t";
	}
	cout << endl;
}

void print_eqlist(PALP::EqList* eqlist) {
	cout << eqlist->ne << endl;
	for (int i = 0; i < eqlist->ne; i++) {
		PALP::Equation eq = eqlist->e[i];
		cout << eq.c << endl;
		for (long w : eq.a) {
			cout << w << "\t";
		}
		cout << endl;
	}

}


int main(int argc, char** argv) {

	// For testing
	/*
	string vertstr1 = "{{1,0,0,0},{1,2,0,0},{0,0,1,0},{2,0,3,4},{-2,0,-1,-2},{-2,-2,-1,-2},{-1,-1,1,0},{-2,-3,-3,-3},{-2,2,-3,-3}}"; //ref
	string vertstr2 = "{{2,0,0,0},{1,5,0,0},{1,0,5,0},{1,0,0,5},{1,1,0,0},{-4,-5,-5,-5}}"; //not
	string vertstr3 = "{{2,0,0,0},{1,6,0,0},{1,0,5,0},{1,0,0,5},{1,3,2,0},{-4,-5,-5,-5}}"; //not
	string vertstr4 = "{{4,-1,-1,-1},{-1,1,0,0},{-1,0,1,0},{-1,0,0,0}}";
	vector<vector<int>> verts = string_to_vertex_list(vertstr1);
	*/
	
	//int polyid = atoi(argv[1]);
	//vector<vector<int>> verts = string_to_vertex_list(argv[2]);

	string line;
	Document doc;
	int polyid;
	vector<vector<int>> verts;
	while (getline(cin, line)) {
		if (doc.Parse(line.c_str()).HasParseError())
	        return 1;
	    assert(doc["POLYID"].IsInt());
	    assert(doc["NVERTS"].IsString());

		polyid = doc["POLYID"].GetInt();
		verts = string_to_vertex_list(doc["NVERTS"].GetString());

		PolyMatrix pmat(verts[0].size(),verts.size());
		for (int i=0; i < verts.size(); i++) {
			for (int j=0; j < verts[0].size(); j++) {
				pmat(j,i) = (double) verts[i][j];
			}
		}

		PolyMatrix origin(verts[0].size(),1);
			for (int i=0; i < verts[0].size(); i++) {
				origin(i,0) = (double) 0;
			}
		
		LatticePolytope p = LatticePolytope(pmat);

		/*cout << "Equations 1:" << endl;
		p.print_equations();*/

		/*
		cout << "The LatticePolytope p has " << p.nvertices() << " vertices" << endl;
		cout << "The LatticePolytope p has ambient dimension " << p.ambient_dim() << endl;
		cout << "The LatticePolytope p has dimension " << p.dim() << endl;
		cout << "The LatticePolytope is ";
		if (!p.is_reflexive()) {
			cout << "not ";
		}
		cout << "reflexive" << endl;
		*/

		/*cout << "Equations 2:" << endl;
		p.print_equations();*/

		/*
		cout << "The polytope has " << p.npoints() << " points" << endl;

		cout << "p.vnlist is: \t";
		PALP::VertexNumList vnlist = p.get_vnlist();
		for (int i = 0; i < p.nvertices(); i++) {
			cout << vnlist.v[i] << "\t";
		}
		cout << endl;

		cout << "The f-vector of p is:" << endl;
		for (int w : p.f_vector()) {
			cout << w << "\t";
		}
		cout << endl;
		*/

		/*cout << "Equations 3:" << endl;
		p.print_equations();*/

		//cout << "The vertex string is: " << p.vertex_string() << endl;

		/*cout << "Equations 4:" << endl;
		p.print_equations();*/

		/*
		if (p.is_reflexive()) {

			cout << "The normal form is:" << endl;
			cout << p.normal_form() << endl;
			cout << "or in string form: " << polymat_to_string(p.normal_form()) << endl;
		

			cout << "The Hodge numbers are:" << endl;
			vector<int> hnums = p.hodge();
			cout << "h11: " << hnums[0] << endl;
			cout << "h21: " << hnums[1] << endl;
		}
		*/

		/*cout << "Equations 5:" << endl;
		p.print_equations();*/

		/*
		PolyMatrix nf = p.normal_form();
		cout << "The normal form is:" << endl;
		cout << nf << endl;
		cout << "or in string form: " << polymat_to_string(nf) << endl;

		cout << "The boundary points are:" << endl;
		cout << p.boundary_points() << endl;
		cout << "of which there are:" << endl;
		cout << p.boundary_points().cols() << endl;

		cout << "The interior points are:" << endl;
		cout << p.interior_points() << endl;
		cout << "of which there are:" << endl;
		cout << p.interior_points().cols() << endl;

		vector<int> orig = {0,0,0,0};
		cout << "p ";
		if (p.interior_contains(orig)) {cout << "contains ";}
		else {cout << "does not contain ";}
		cout << "the origin" << endl;
		*/

		/*vector<int> pt1 = {1,4,1,0};
		vector<int> pt2 = {0,-1,-1,0};
		vector<int> pt3 = {2,4,6,8};
		cout << p.boundary_contains(pt1) << endl;
		cout << p.boundary_contains(pt2) << endl;
		cout << p.boundary_contains(pt3) << endl;
		cout << p.interior_contains(pt1) << endl;
		cout << p.interior_contains(pt2) << endl;
		cout << p.interior_contains(pt3) << endl;*/

		if (p.is_reflexive()) {
			// For the dual
			LatticePolytope dp = p.polar();

			/*
			cout << "The dual LatticePolytope dp has " << dp.nvertices() << " vertices" << endl;
			cout << "The dual LatticePolytope dp has ambient dimension " << dp.ambient_dim() << endl;
			cout << "The dual LatticePolytope dp has dimension " << dp.dim() << endl;
			cout << "The dual LatticePolytope is ";
			if (!dp.is_reflexive()) {
				cout << "not ";
			}
			cout << "reflexive" << endl;

			cout << "The dual polytope has " << dp.npoints() << " points" << endl;

			cout << "dp.vnlist is: \t";
			PALP::VertexNumList vnlistd = dp.get_vnlist();
			for (int i = 0; i < dp.nvertices(); i++) {
				cout << vnlistd.v[i] << "\t";
			}
			cout << endl;

			if (dp.is_reflexive()) {
				cout << "The dual normal form is:" << endl;
				cout << dp.normal_form() << endl;
				cout << "or in string form: " << polymat_to_string(dp.normal_form()) << endl;
			

				cout << "The dual Hodge numbers are:" << endl;
				vector<int> dhnums = dp.hodge();
				cout << "h11: " << dhnums[0] << endl;
				cout << "h21: " << dhnums[1] << endl;
			}
			*/

			//cout << "The vertices of dp are:" << endl;
			cout << "set POLY {\"POLYID\":" << polyid << "} {\"DVERTS\":\"" << polymat_to_string(dp.vertices()) << "\",";

			//cout << "The points of the 2-skeleton of dp are:" << endl;
			cout << "\"DRESVERTS\":\"" << polymat_to_string(dp.p_boundary_points(3)) << "\",";

			/*
			cout << "The integral points of dp are: " << endl;
			cout << dp.integral_points() << endl;

			cout << "The 3d faces of dp are: " << endl;
			*/

			vector<PolyMatrix> dfaces3 = dp.faces(3);
			
			/*
			for (PolyMatrix fm : dfaces3) {
				cout << "=====" << endl;
				cout << fm << endl;
			}

			cout << "Boundary and interior points of faces of dp are: " << endl;

			for (PolyMatrix fm : dfaces3) {
				LatticePolytope fmp = LatticePolytope(add_col(fm, origin));
				PolyMatrix fmverts = remove_col(fmp.vertices(), origin);
				PolyMatrix fm1bdry = remove_col(fmp.p_boundary_points(1), origin);
				PolyMatrix fm1inter = remove_col(fmp.p_interior_points(1), origin);
				PolyMatrix fm2bdry = remove_col(fmp.p_boundary_points(2), origin);
				PolyMatrix fm2inter = remove_col(fmp.p_interior_points(2), origin);
				PolyMatrix fm3bdry = remove_col(fmp.p_boundary_points(3), origin);
				PolyMatrix fm3inter = remove_col(fmp.p_interior_points(3), origin);
				PolyMatrix fmpoints = remove_col(fmp.integral_points(), origin);
				
				cout << "=====" << endl;
				cout << fmverts.cols() << " " << fm1bdry.cols() << " " << fm1inter.cols() << " "  << fm2bdry.cols() << " " << fm2inter.cols()  << " "  << fm3bdry.cols() << " " << fm3inter.cols() << " " << fmpoints.cols() << endl;
			}

			cout << "The 2-skeletons of the 3d faces of dp are: " << endl;
			*/

			vector<PolyMatrix> nffs;
			vector<PolyMatrix> nf2skels;
			for (PolyMatrix fm : dfaces3) {
				LatticePolytope fmp = LatticePolytope(add_col(fm, origin));

				PolyMatrix nf = fmp.normal_form();
				PolyMatrix nff = remove_col(nf, origin);
				nffs.push_back(nff);

				LatticePolytope nfp = LatticePolytope(nf);

				PolyMatrix nf2skel = remove_col(nfp.p_boundary_points(3), origin);
				nf2skels.push_back(nf2skel);
				//cout << "=====" << endl;
				//cout << nf2skel << endl;
			}

			vector<PolyMatrix> nffsuniq;
			vector<PolyMatrix> nf2skelsuniq;
			vector<int> nffsmult;
			int uniq;
			int equal;
			for (int i = 0; i < nffs.size(); i++) {
				uniq = 1;
				for (int j = 0; j < nffsuniq.size(); j++) {
					if (nffs[i].rows() == nffsuniq[j].rows() && nffs[i].cols() == nffsuniq[j].cols()) {
						equal = 1;
						for (int k = 0; k < nffs[i].rows(); k++) {
							for (int l = 0; l < nffs[i].cols(); l++) {
								if (nffs[i](k,l) != nffsuniq[j](k,l)) {
									equal = 0;
								}
							}
						}
						if (equal == 1) {
							uniq = 0;
							nffsmult[j]++;
						}
					}
				}
				if (uniq == 1) {
					nffsuniq.push_back(nffs[i]);
					nf2skelsuniq.push_back(nf2skels[i]);
					nffsmult.push_back(1);
				}
			}

			//cout << "The unique 2-skeletons of the 3d faces of dp are: " << endl;
			cout << "\"FACETLIST\":[";

			for (int i = 0; i < nffsuniq.size(); i++) {
				//cout << "=====" << endl;
				if (i > 0) {
					cout << ",";
				}
				cout << "{\"NFORM\":\"" << polymat_to_string(nffsuniq[i]) << "\",\"NINST\":" << nffsmult[i] << "}";
			}
			cout << "]}";

			cout << endl;

			for (int i = 0; i < nffsuniq.size(); i++) {
				//cout << "=====" << endl;
				cout << "set FACET {\"NFORM\":\"" << polymat_to_string(nffsuniq[i]) << "\"} {\"NFORM2SKEL\":\"" << polymat_to_string(nf2skelsuniq[i]) << "\",\"facetfinetriangsMARK\":false,\"facetallfinetriangsMARK\":false}" << endl;
			}

			cout << endl;

			/*
			cout << "The f-vector of dp is:" << endl;
			for (int w : dp.f_vector()) {
				cout << w << "\t";
			}
			cout << endl;

			cout << "The dual vertex string is: " << dp.vertex_string() << endl;

			cout << "The boundary points are:" << endl;
			cout << dp.boundary_points() << endl;
			cout << "of which there are:" << endl;
			cout << dp.boundary_points().cols() << endl;

			cout << "The interior points are:" << endl;
			cout << dp.interior_points() << endl;
			cout << "of which there are:" << endl;
			cout << dp.interior_points().cols() << endl;

			cout << "dp ";
			if (dp.interior_contains(orig)) {cout << "contains ";}
			else {cout << "does not contain ";}
			cout << "the origin" << endl;	
			*/																
		}
	}
}