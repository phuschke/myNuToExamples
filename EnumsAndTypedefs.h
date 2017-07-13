#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Core>
#include "mechanics/MechanicsEnums.h"
#include "visualize/VisualizeEnum.h"


using std::cout;
using std::endl;
using NuTo::Constitutive::ePhaseFieldEnergyDecomposition;
using NuTo::Constitutive::eConstitutiveType;
using NuTo::Constitutive::eConstitutiveParameter;
using NuTo::Node::eDof;
using NuTo::Interpolation::eTypeOrder;
using NuTo::Interpolation::eShapeType;
using NuTo::eGroupId;
using NuTo::eVisualizeWhat;

using NuTo::eDirection;

using Eigen::VectorXd;
using Eigen::MatrixXd;

using Eigen::Vector2d;
using Eigen::Matrix2d;