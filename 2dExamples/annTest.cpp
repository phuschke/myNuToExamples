#include "nuto/math/FullMatrix.h"
#include "nuto/mechanics/structures/unstructured/Structure.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStrain3D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStress3D.h"
#include "nuto/mechanics/constitutive/mechanics/GradientDamageEngineeringStress.h"
#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/constitutive/ConstitutiveStaticDataBase.h"
#include "nuto/mechanics/constitutive/mechanics/ConstitutiveStaticDataGradientDamage1D.h"

#include "nuto/math/SparseMatrixCSRVector2General.h"

#include "nuto/mechanics/timeIntegration/NewmarkDirect.h"
#include "nuto/base/Logger.h"

#include <boost-1_55/boost/filesystem.hpp>
#include <boost-1_55/boost/lexical_cast.hpp>
#include <boost-1_55/boost/progress.hpp>

#include <iostream>
#include <fstream>
#include <string>
#include <array>
#include <cmath>

#include "ANN/ANN.h"






void printPt(ANNpoint p)          // print point
{
//    int dim = 3; // dimension
//    std::cout << "(" << p[0];
//    for (int i = 1; i < dim; i++) {
//        std::cout << ", " << p[i];
//    }
//    std::cout << ")\n";
}




int main()
{
    int k = 3; // number of nearest neighbors
    int dim = 3; // dimension
    double eps = 1e-6; // error bound
    int maxPts = 10; // maximum number of data points
int numPoints;







//ANNcoord* queryPoint;
//ANNpoint* dataPoints;
ANNdist* dists;


ANNidx* nnIds;

ANNkd_tree* kdTree;

ANNcoord* queryPoint = annAllocPt(dim);
ANNpoint* dataPoints = annAllocPts(maxPts, dim);

queryPoint[0] = 1.0;
queryPoint[1] = 1.0;
queryPoint[2] = 1.0;



Eigen::VectorXd bla(3);
bla(0,0) = 6;
bla(1,0) = 1;
bla(2,0) = 3;

dataPoints[1] = const_cast<double*>(bla.data());


nnIds = new ANNidx[k];
dists = new ANNdist[k];

numPoints = 3;

kdTree = new ANNkd_tree(dataPoints, numPoints, dim);


kdTree->annkSearch(queryPoint, k,nnIds, dists, eps);

std::cout << "NN: Index Distance\n";
for (int i = 0; i < k; i++) { // print summary
dists[i] = sqrt(dists[i]); // unsquare distance
std::cout << i << " " << nnIds[i] << " " << dists[i] << "\n";
}


delete [] nnIds; // clean things up
delete [] dists;
delete kdTree;

delete [] queryPoint;
delete [] dataPoints[0];
delete [] dataPoints;


annClose(); // done with ANN
//delete queryPoint;
//delete [] dataPoints;



//Eigen::VectorXd bla(3);
//bla(0,0) = 1;
//bla(1,0) = 2;
//bla(2,0) = 3;
//
//
//
//double* p = new double[3];
//
//p = bla.data();
//
//for ( int i = 0; i < 5; i++ )
//   {
//       std::cout << "*(p + " << i << ") : ";
//       std::cout << *(p + i) << std::endl;
//   }


}
