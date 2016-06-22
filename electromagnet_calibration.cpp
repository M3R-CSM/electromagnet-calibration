#include "electromagnet_calibration.h"

#include <iomanip>
#include <Eigen/Jacobi>

ElectromagnetCalibration::MagneticWorkSpace::MagneticWorkSpace()
{
    xMin = xMax = yMin = yMax = zMin = zMax = 0;
}

ElectromagnetCalibration::MagneticWorkSpace::MagneticWorkSpace(double size)
{
    xMin = -size;
    xMax = size;
    yMin = -size;
    yMax = size;
    zMin = -size;
    zMax = size;
}


ElectromagnetCalibration::MagneticWorkSpace::MagneticWorkSpace(double xMin_, double xMax_, double yMin_, double yMax_, double zMin_, double zMax_ )
{
    xMin = xMin_;
    xMax = xMax_;
    yMin = yMin_;
    yMax = yMax_;
    zMin = zMin_;
    zMax = zMax_;
}

MagneticMeasurement::MagneticMeasurement()
{
    Field.setZero(3,1);
    Position.setZero(3,1);
    AppliedCurrentVector.setZero(0,1);

}

MagneticMeasurement::MagneticMeasurement(const Eigen::Vector3d &F, const Eigen::Vector3d &P, const Eigen::VectorXd &C)
{
    Field = F;
    Position = P;
    AppliedCurrentVector = C;
}


ElectromagnetCalibration::ElectromagnetCalibration( std::string calibrationFileName )
{
    bool calibrationFileLoadSucessful = loadCalibration(calibrationFileName);

    assert( calibrationFileLoadSucessful );

    checkSourcePositions();
}

ElectromagnetCalibration::ElectromagnetCalibration()
{
    /*< Default Constructor only available to inheriting classes **/
    this->name = "NONE";
    this->coilList.clear();
    this->use_offset = false;

    checkSourcePositions();
}

//
// \brief the constructor from a coil list.
// \param coilList The list of coils and their respective sources.
// \param dc_field_offset The dc offset field, if any.
//
ElectromagnetCalibration::ElectromagnetCalibration(std::string systemName_, const MagneticWorkSpace& workSpace_, const std::vector<ScalorPotential>& coilList_, const ScalorPotential& dc_field_offset_ )
{
    /*< Default Constructor only available to inheriting classes **/
    this->coilList = coilList_;
    this->name = systemName_ ;
    this->offset = dc_field_offset_;
    this->use_offset = ( 0 != offset.getNumberOfSources());
    this->workSpace = workSpace_;

    checkSourcePositions();
}

// The following functions return the field and/or gradient given a current vector and a position.
//   If no position is provided, it is assumed to be at the workspace origin. For the combined vector,
//   the gradient matrix has been repacked into a 5 element vetor form (because it is symetric and has zero trace).
//   The order of the gradient terms is: [dBx/dx, dBx/dy, dBx/dz, dBy/dy, dBy/dz].
Eigen::Vector3d ElectromagnetCalibration::fieldAtPoint(    const Eigen::VectorXd& currentVector, const Eigen::Vector3d& position ) const
{
    assert( currentVector.size() == coilList.size() );
    Eigen::Vector3d field(0,0,0);
    std::vector<ScalorPotential>::const_iterator coilIT = coilList.begin();
    int coilNum = 0;
    for( ; coilIT != coilList.end(); coilIT++, coilNum ++ )
    {
        field += coilIT->getGradient(position)*currentVector(coilNum);
    }
    if( use_offset )
    {
        field += offset.getGradient(position);
    }


    return field;
}

Eigen::Matrix3d ElectromagnetCalibration::gradientAtPoint( const Eigen::VectorXd& currentVector, const Eigen::Vector3d& position ) const
{

    assert( currentVector.size() == coilList.size() ); // make sure there are as many currents as coils
    Eigen::MatrixXd actuationMatrix = gradientCurrentJacobian( position );

    return remapGradientVector(actuationMatrix*currentVector + offsetFieldAndGradientAtPoint(position).tail<5>());

}

Vector8d ElectromagnetCalibration::fieldAndGradientAtPoint( const Eigen::VectorXd& currentVector, const Eigen::Vector3d& position ) const
{
    assert( currentVector.size() == coilList.size() ); // make sure there are as many currents as coils
    Eigen::MatrixXd actuationMatrix = fieldAndGradientCurrentJacobian( position );

    return actuationMatrix*currentVector + offsetFieldAndGradientAtPoint(position);
}

Vector8d ElectromagnetCalibration::offsetFieldAndGradientAtPoint( const Eigen::Vector3d& position ) const
{
    Vector8d offsetFieldAndGradient = Vector8d::Zero();

    if( use_offset )
    {

        ScalorPotentialState fieldInfo = offset.getState(position);


        offsetFieldAndGradient.head<3>() = fieldInfo.firstSpatialDerivative;
        offsetFieldAndGradient.tail<5>() = remapGradientMatrix(fieldInfo.secondSpatialDerivative);
    }

    return offsetFieldAndGradient;
}


// The following functions return the Current Jacobian of the field and/or gradient.  This matrix can be inverted to determine the
//   currents necessary to achieve a desired field and gradient.  The gradient portion of the matrix has been vectorized into 5 elements.
//   They are [dBx/dx, dBx/dy, dBx/dz, dBy/dy, dBy/dz].
Eigen::MatrixXd ElectromagnetCalibration::fieldCurrentJacobian( const Eigen::Vector3d& position ) const
{
    Eigen::MatrixXd actuationMatrix = Eigen::MatrixXd::Zero(3,coilList.size());



    if( !pointInWorkspace(position) )
    {
        cout << "Warning: Requesting point out side of magnetic workspace. " << position.transpose() << " (Xmin,Xmax): (" << this->workSpace.xMin << ", " << this->workSpace.xMax << ") " <<
                " (Ymin,Ymax): (" << this->workSpace.yMin << ", " << this->workSpace.yMax << ") " <<
                " (Zmin,Zmax): (" << this->workSpace.zMin << ", " << this->workSpace.zMax << ") " << endl;
    }

    std::vector<ScalorPotential>::const_iterator coilIT = coilList.begin();
    int coilNum = 0;
    for( ; coilIT != coilList.end(); coilIT++, coilNum ++ )
    {
        ScalorPotentialState fieldInfo = coilIT->getState(position);
        actuationMatrix.block<3,1>(0,coilNum) = fieldInfo.firstSpatialDerivative;

    }

    return actuationMatrix;
}

Eigen::MatrixXd ElectromagnetCalibration::gradientCurrentJacobian( const Eigen::Vector3d& position ) const
{

    Eigen::MatrixXd actuationMatrix = Eigen::MatrixXd::Zero(5,coilList.size());


    if( !pointInWorkspace(position) )
    {
        cout << "Warning: Requesting point out side of magnetic workspace. " << position.transpose() << " (Xmin,Xmax): (" << this->workSpace.xMin << ", " << this->workSpace.xMax << ") " <<
                " (Ymin,Ymax): (" << this->workSpace.yMin << ", " << this->workSpace.yMax << ") " <<
                " (Zmin,Zmax): (" << this->workSpace.zMin << ", " << this->workSpace.zMax << ") " << endl;
    }

    std::vector<ScalorPotential>::const_iterator coilIT = coilList.begin();
    int coilNum = 0;
    for( ; coilIT != coilList.end(); coilIT++, coilNum ++ )
    {
        ScalorPotentialState fieldInfo = coilIT->getState(position);
        actuationMatrix.block<5,1>(0,coilNum) = remapGradientMatrix(fieldInfo.secondSpatialDerivative);

    }

    return actuationMatrix;
}

Eigen::MatrixXd ElectromagnetCalibration::fieldAndGradientCurrentJacobian( const Eigen::Vector3d& position ) const
{

    Eigen::MatrixXd actuationMatrix = Eigen::MatrixXd::Zero(8,coilList.size());


    if( !pointInWorkspace(position) )
    {
        cout << "Warning: Requesting point out side of magnetic workspace. " << position.transpose() << " (Xmin,Xmax): (" << this->workSpace.xMin << ", " << this->workSpace.xMax << ") " <<
                " (Ymin,Ymax): (" << this->workSpace.yMin << ", " << this->workSpace.yMax << ") " <<
                " (Zmin,Zmax): (" << this->workSpace.zMin << ", " << this->workSpace.zMax << ") " << endl;
    }

    std::vector<ScalorPotential>::const_iterator coilIT = coilList.begin();
    int coilNum = 0;
    for( ; coilIT != coilList.end(); coilIT++, coilNum ++ )
    {
        ScalorPotentialState fieldInfo = coilIT->getState(position);
        actuationMatrix.block<3,1>(0,coilNum) = fieldInfo.firstSpatialDerivative;
        actuationMatrix.block<5,1>(3,coilNum) = remapGradientMatrix(fieldInfo.secondSpatialDerivative);

    }

    return actuationMatrix;
}

// This function returns the 5x3 Jacobian describing how gradient changes with position.  The gradient change is a 5x3 packing of a 3x3x3 tensor.
//  The first column is how the gradient vector packing [dBx/dx, dBx/dy, dBx/dz, dBy/dy, dBy/dz] changes with x, the second is how it changes with y, and third is
//  how it changes with z.
Eigen::Matrix<double,5,3> ElectromagnetCalibration::gradientPositionJacobian( const Eigen::VectorXd& currentVector, const Eigen::Vector3d& position ) const
{

    assert( currentVector.size() == coilList.size() );// make sure there are as many currents as coils

    Eigen::Matrix<double,5,3> gradJacobian = Eigen::MatrixXd::Zero(5,3);

    if( currentVector.norm() == 0 && !useOffset() )
    {
        // NOthing to do!
        return gradJacobian;
    }

    if( !pointInWorkspace(position) )
    {
        cout << "Warning: Requesting point out side of magnetic workspace. " << position.transpose() << " (Xmin,Xmax): (" << this->workSpace.xMin << ", " << this->workSpace.xMax << ") " <<
                " (Ymin,Ymax): (" << this->workSpace.yMin << ", " << this->workSpace.yMax << ") " <<
                " (Zmin,Zmax): (" << this->workSpace.zMin << ", " << this->workSpace.zMax << ") " << endl;
    }

    std::vector<ScalorPotential>::const_iterator coilIT = coilList.begin();
    int coilNum = 0;
    for( ; coilIT != coilList.end(); coilIT++, coilNum ++ )
    {
        ScalorPotentialState fieldInfo = coilIT->getState(position);
        gradJacobian += fieldInfo.thirdSpatialDerivative*currentVector(coilNum);

    }

    if( use_offset )
    {
        ScalorPotentialState fieldInfo = offset.getState(position);
        gradJacobian += fieldInfo.thirdSpatialDerivative;
    }

    return gradJacobian;
}


// This function returns the field, gradient, gradientJacobian, and field/gradient current jacobain.  It is more efficient than requesting them seporately
void ElectromagnetCalibration::fullMagneticState( Eigen::Vector3d& fieldAtPoint, Eigen::Matrix<double,8,3>& fieldGradientPositionJacobian, Eigen::MatrixXd& fieldGradientCurrentJacobian, const Eigen::VectorXd& currentVector, const Eigen::Vector3d& position ) const
{
    assert( currentVector.size() == coilList.size() );// make sure there are as many currents as coils

    fieldAtPoint.setZero();
    fieldGradientPositionJacobian.setZero();
    fieldGradientCurrentJacobian.setZero(8,coilList.size());

    Vector5d gradAtPt;
    gradAtPt.setZero();


    if( currentVector.norm() == 0 && !this->useOffset() )
    {
        // NOthing to do!
        return ;
    }

    if( !pointInWorkspace(position) )
    {
        cout << "Warning: Requesting point out side of magnetic workspace. " << position.transpose() << " (Xmin,Xmax): (" << this->workSpace.xMin << ", " << this->workSpace.xMax << ") " <<
                " (Ymin,Ymax): (" << this->workSpace.yMin << ", " << this->workSpace.yMax << ") " <<
                " (Zmin,Zmax): (" << this->workSpace.zMin << ", " << this->workSpace.zMax << ") " << endl;
    }

    std::vector<ScalorPotential>::const_iterator coilIT = coilList.begin();
    int coilNum = 0;
    for( ; coilIT != coilList.end(); coilIT++, coilNum ++ )
    {
        ScalorPotentialState fieldInfo = coilIT->getState(position);

        fieldGradientCurrentJacobian.block<3,1>(0,coilNum) = fieldInfo.firstSpatialDerivative;
        fieldGradientCurrentJacobian.block<5,1>(3,coilNum) = remapGradientMatrix(fieldInfo.secondSpatialDerivative);

        fieldGradientPositionJacobian.block<5,3>(3,0) += fieldInfo.thirdSpatialDerivative*currentVector(coilNum);
        fieldGradientPositionJacobian.block<3,3>(0,0) += fieldInfo.secondSpatialDerivative*currentVector(coilNum);
        fieldAtPoint += fieldInfo.firstSpatialDerivative*currentVector(coilNum);
    }

    if( use_offset )
    {
        ScalorPotentialState fieldInfo = offset.getState(position);

        fieldGradientPositionJacobian.block<5,3>(3,0) += fieldInfo.thirdSpatialDerivative;
        fieldGradientPositionJacobian.block<3,3>(0,0) += fieldInfo.secondSpatialDerivative;
        fieldAtPoint += fieldInfo.firstSpatialDerivative;

    }

}

MagneticState ElectromagnetCalibration::fullMagneticState(const Eigen::VectorXd& currentVector, const Eigen::Vector3d& position ) const
{
    assert( currentVector.size() == coilList.size() );// make sure there are as many currents as coils

    MagneticState returnState;
    returnState.FieldGradientActuationMatrix.setZero(8,coilList.size());


    if( !pointInWorkspace(position) )
    {
        cout << "Warning: Requesting point out side of magnetic workspace. " << position.transpose() << " (Xmin,Xmax): (" << this->workSpace.xMin << ", " << this->workSpace.xMax << ") " <<
                " (Ymin,Ymax): (" << this->workSpace.yMin << ", " << this->workSpace.yMax << ") " <<
                " (Zmin,Zmax): (" << this->workSpace.zMin << ", " << this->workSpace.zMax << ") " << endl;
    }

    std::vector<ScalorPotential>::const_iterator coilIT;
    int coilNum;
    for( coilIT = coilList.begin(), coilNum = 0;
         coilIT != coilList.end();
         coilIT++, coilNum ++ )
    {
        ScalorPotentialState fieldInfo = coilIT->getState(position);

        returnState.FieldGradientActuationMatrix.block<3,1>(0,coilNum) = fieldInfo.firstSpatialDerivative;
        returnState.FieldGradientActuationMatrix.block<5,1>(3,coilNum) =  remapGradientMatrix(fieldInfo.secondSpatialDerivative);

        returnState.Field += fieldInfo.firstSpatialDerivative * currentVector(coilNum);
        returnState.GradientPositionJacobian += fieldInfo.thirdSpatialDerivative*currentVector(coilNum);
        returnState.Gradient += fieldInfo.secondSpatialDerivative*currentVector(coilNum);

    }

    if( use_offset )
    {
        ScalorPotentialState fieldInfo = offset.getState(position);

        returnState.Field += fieldInfo.firstSpatialDerivative;
        returnState.GradientPositionJacobian += fieldInfo.thirdSpatialDerivative;
        returnState.Gradient += fieldInfo.secondSpatialDerivative;
    }

    return returnState;
}



// These two functions convert a vector packing of the gradient to a matrix packing and vice versa.
Eigen::Matrix3d ElectromagnetCalibration::remapGradientVector(const Vector5d& gradVector )
{
    Eigen::Matrix3d gradientMatrix;
    // repack vector gradient into matrix gradient
    gradientMatrix.leftCols<1>() = gradVector.topRows<3>();
    gradientMatrix.block<2,1>(1,1) = gradVector.bottomRows<2>();
    gradientMatrix(0,1) = gradientMatrix(1,0);
    gradientMatrix(0,2) = gradientMatrix(2,0);
    gradientMatrix(1,2) = gradientMatrix(2,1);
    gradientMatrix(2,2) = -1.0*(gradientMatrix(0,0)+gradientMatrix(1,1));

    return gradientMatrix;
}

Vector5d ElectromagnetCalibration::remapGradientMatrix( const Eigen::Matrix3d& gradMatrix)
{
    Vector5d gradVec;
    gradVec.topRows<3>() = gradMatrix.leftCols<1>();
    gradVec(3) = gradMatrix(1,1);
    gradVec(4) = gradMatrix(2,1);
    return gradVec;
}

Eigen::MatrixXd ElectromagnetCalibration::packForceMatrix(const Eigen::Vector3d& moment)
{
    Eigen::MatrixXd ForceMatrix(3,5);
    ForceMatrix <<  moment(0),  moment(1), moment(2),  0,         0,
            0,          moment(0), 0,          moment(1), moment(2),
            -moment(2), 0,         moment(0), -moment(2), moment(1);
    return ForceMatrix;
}

bool ElectromagnetCalibration::loadCalibration(std::string fileName)
{
    coilList.clear();

    // reinitialize offset
    offset = ScalorPotential();

    YAML::Node systemDefinition = YAML::LoadFile(fileName);
    this->name = systemDefinition["System_Name"].as<std::string>();
    int numCoils = systemDefinition["Coil_List"].size();
    this->workSpace.xMin = systemDefinition["Workspace_Dimensions"][0][0].as<double>();
    this->workSpace.xMax = systemDefinition["Workspace_Dimensions"][0][1].as<double>();
    this->workSpace.yMin = systemDefinition["Workspace_Dimensions"][1][0].as<double>();
    this->workSpace.yMax = systemDefinition["Workspace_Dimensions"][1][1].as<double>();
    this->workSpace.zMin = systemDefinition["Workspace_Dimensions"][2][0].as<double>();
    this->workSpace.zMax = systemDefinition["Workspace_Dimensions"][2][1].as<double>();

    for( unsigned int coil = 0; coil<numCoils; coil++ )
    {
        std::string coilTag = systemDefinition["Coil_List"][coil].as<std::string>();

        ScalorPotential newCoil;

        int numSources = systemDefinition[coilTag]["Source_List"].size();

        for( unsigned int src = 0; src < numSources; src ++ )
        {
            ScalorPotential::srcStruct newSrc;
            std::string srcTag = systemDefinition[coilTag]["Source_List"][src].as<std::string>();

            std::vector<double> coeff = systemDefinition[coilTag][srcTag]["A_Coeff"].as< std::vector<double> >();
            for( unsigned int i=0; i<coeff.size(); i++ )
                newSrc.A_Coeff.push_back(ScalorPotential::srcCoeff(coeff[i],i+1));

            coeff = systemDefinition[coilTag][srcTag]["B_Coeff"].as< std::vector<double> >();
            for( unsigned int i=0; i<coeff.size(); i++ )
                newSrc.B_Coeff.push_back(ScalorPotential::srcCoeff(coeff[i],i+1));

            newSrc.srcDirection = systemDefinition[coilTag][srcTag]["Source_Direction"].as< Eigen::Vector3d >();
            newSrc.srcDirection.normalize(); // normalize in place
            newSrc.srcPosition = systemDefinition[coilTag][srcTag]["Source_Position"].as< Eigen::Vector3d >();

            newCoil.setSourceStruct(src,newSrc); // add to coil's source list
        }

        if( coilTag.find("Offset") == std::string::npos )
        {
            // "Offset" not found
            this->coilList.push_back(newCoil); // add to coil List
        } else
        {
            // Offset tag is found
            offset = newCoil;
        }
    }

    use_offset = offset.getNumberOfSources();
    return coilList.size() > 0 || use_offset;
}

///@brief  writes a new calibration file
/// @param a string pointing to the location for the yaml formated calbiration file
bool ElectromagnetCalibration::writeCalibration(std::string fileName) const
{
    YAML::Node systemDefinition;
    systemDefinition["System_Name"] = this->name;

    std::vector<ScalorPotential>::const_iterator coilIT = coilList.begin();
    unsigned int coilNum =0;
    for( ; coilIT != coilList.end(); coilIT ++, coilNum ++ )
    {
        stringstream coilName;
        coilName << "Coil_" << coilNum;
        systemDefinition["Coil_List"].push_back(coilName.str());

        for( unsigned int srcNum = 0; srcNum < coilIT->getNumberOfSources(); srcNum ++ )
        {
            ScalorPotential::srcStruct src = coilIT->getSourceStruct(srcNum);
            stringstream srcName;
            srcName << "Src_" << srcNum;
            systemDefinition[coilName.str()]["Source_List"].push_back(srcName.str());
            std::vector<double> coeff;
            for( unsigned int i=0; i<src.A_Coeff.size(); i++ )
                coeff.push_back(src.A_Coeff[i].coeff);
            systemDefinition[coilName.str()][srcName.str()]["A_Coeff"] = coeff;

            coeff.clear();
            for( unsigned int i=0; i<src.B_Coeff.size(); i++ )
                coeff.push_back(src.B_Coeff[i].coeff);
            systemDefinition[coilName.str()][srcName.str()]["B_Coeff"] = coeff;
            systemDefinition[coilName.str()][srcName.str()]["Source_Direction"] = src.srcDirection;
            systemDefinition[coilName.str()][srcName.str()]["Source_Position"] = src.srcPosition;
        }
    }

    if( this->hasOffset() )
    {
        stringstream coilName;
        coilName << "Offset";
        systemDefinition["Coil_List"].push_back(coilName.str());

        for( unsigned int srcNum = 0; srcNum < offset.getNumberOfSources(); srcNum ++ )
        {
            ScalorPotential::srcStruct src = offset.getSourceStruct(coilNum);
            stringstream srcName;
            srcName << "Src_" << srcNum;
            systemDefinition[coilName.str()]["Source_List"].push_back(srcName.str());
            std::vector<double> coeff;
            for( unsigned int i=0; i<src.A_Coeff.size(); i++ )
                coeff.push_back(src.A_Coeff[i].coeff);
            systemDefinition[coilName.str()][srcName.str()]["A_Coeff"] = coeff;

            coeff.clear();
            for( unsigned int i=0; i<src.A_Coeff.size(); i++ )
                coeff.push_back(src.B_Coeff[i].coeff);
            systemDefinition[coilName.str()][srcName.str()]["B_Coeff"] = coeff;
            systemDefinition[coilName.str()][srcName.str()]["Source_Direction"] = src.srcDirection;
            systemDefinition[coilName.str()][srcName.str()]["Source_Position"] = src.srcPosition;
        }
    }

    systemDefinition["Workspace_Dimensions"][0][0] = this->workSpace.xMin;
    systemDefinition["Workspace_Dimensions"][0][1] = this->workSpace.xMax;
    systemDefinition["Workspace_Dimensions"][1][0] = this->workSpace.yMin;
    systemDefinition["Workspace_Dimensions"][1][1] = this->workSpace.yMax;
    systemDefinition["Workspace_Dimensions"][2][0] = this->workSpace.zMin;
    systemDefinition["Workspace_Dimensions"][2][1] = this->workSpace.zMax;


    std::ofstream fout;
    fout.open(fileName.c_str(),std::ofstream::out|std::ofstream::trunc);
    fout << systemDefinition;
    fout.close();

    return fout.good();
}

int ElectromagnetCalibration::getNumberOfCoils() const
{
    return coilList.size();
}

///@brief returns the number of sources for the given coil
int ElectromagnetCalibration::getNumberOfSources( unsigned int coilNum ) const
{
    if( coilNum >= (coilList.size() + useOffset()) )
        return 0;
    else if( coilNum == coilList.size() )
        return offset.getNumberOfSources();
    else
        return coilList[coilNum].getNumberOfSources();
}

///@brief returns the number of coefficients for the given source
int ElectromagnetCalibration::getNumberOfCoeffients( unsigned int coilNum, unsigned int srcNum ) const
{
    if( coilNum >= (coilList.size() + useOffset()) )
        return 0;
    else if( coilNum == coilList.size() )
    {
        if( srcNum < offset.getNumberOfSources() )
            return offset.getSourceStruct(srcNum).A_Coeff.size() + offset.getSourceStruct(srcNum).B_Coeff.size();
    }
    else if( srcNum < coilList[coilNum].getNumberOfSources() )
        return coilList[coilNum].getSourceStruct(srcNum).A_Coeff.size() + coilList[coilNum].getSourceStruct(srcNum).B_Coeff.size();

    return 0;

}

bool ElectromagnetCalibration::hasOffset() const
{
    return offset.getNumberOfSources();
}

std::string ElectromagnetCalibration::getName() const
{
    return name;
}

bool ElectromagnetCalibration::pointInWorkspace( const Eigen::Vector3d& position ) const
{
    return workSpace.xMin <= position(0) && workSpace.xMax >= position(0)
            && workSpace.yMin <= position(1) && workSpace.yMax >= position(1)
            && workSpace.zMin <= position(2) && workSpace.zMax >= position(2);
}

void ElectromagnetCalibration::setWorkSpace(const ElectromagnetCalibration::MagneticWorkSpace& ws)
{
    workSpace = ws;
    return;
}

ElectromagnetCalibration::MagneticWorkSpace ElectromagnetCalibration::getWorkSpace() const
{
    return workSpace;
}

void ElectromagnetCalibration::useOffset(bool offsetOn)
{
    use_offset = offsetOn;
}

void ElectromagnetCalibration::useOffset(const ScalorPotential& newOffset )
{
    offset = newOffset;
    use_offset = true;
    return;
}

bool ElectromagnetCalibration::useOffset() const
{
    return use_offset;
}

bool ElectromagnetCalibration::checkSourcePositions(bool printWarning) const
{
    std::vector<ScalorPotential>::const_iterator coilIterator;
    std::vector<ScalorPotential::srcStruct>::const_iterator srcIterator;
    int numTooClose = 0;
    int numTooFar = 0;
    // normalize directions, check positions relative to workspace, and source counts
    for( coilIterator  = coilList.begin();
         coilIterator != coilList.end();
         coilIterator ++)
    {
        for(srcIterator  = coilIterator->srcList.begin();
            srcIterator != coilIterator->srcList.end();
            srcIterator ++ )
        {

            double distanceSq = (srcIterator->srcPosition - pCenter).squaredNorm();

            if( distanceSq < rMinSq)
            {
                numTooClose ++;
            }else if( distanceSq > rMaxSq )
            {
                numTooFar ++;
            }
        }
    }

    if(printWarning && (numTooClose+numTooFar) )
    {
        cout << "There are " << numTooClose << " sources too close to the workspace center and " << numTooFar << " too far from the workspace center." << endl;
    }



    return (numTooClose+numTooFar == 0);
}



void ElectromagnetCalibration::calibrate(std::string calibrationName, const std::vector<MagneticMeasurement>& dataList, bool printProgress, bool printStats_, calibration_constraints constraint,double minimumSourceToCenterDistance, double maximumSourceToCenterDistance, double converganceTolerance, int maxAttempts, int numberOfConvergedIterations)
{
    name = calibrationName;

    if( constraint == HEADING_AND_POSITION )
        minRadIsActive = true;
    else
        minRadIsActive = false;
    nConst = 2;

    posWeight = 1;

    double sqrtEpsilon = std::sqrt(std::numeric_limits<double>::epsilon());
    double rmsError_this, rmsError_last;
    int iteration = 0;
    bool converged;
    int convergedIterations = 0;


    do{
        numberOfMeasurements = dataList.size();
        numberOfParameters = 0;
        numberOfConstraints = 0;
        numberOfSources = 0;

        std::vector<ScalorPotential>::iterator coilIT = coilList.begin();
        std::vector<MagneticMeasurement>::const_iterator measIT;

        // count problem size
        int coilNum = 0;
        for( ; coilIT != coilList.end(); coilIT ++, coilNum ++ )
        {
            numberOfSources += coilIT->getNumberOfSources();
            numberOfConstraints += nConst*coilIT->getNumberOfSources();
            numberOfParameters += coilIT->getNumCalibrationParameters();
        }
        if( useOffset() )
        {
            numberOfSources += offset.getNumberOfSources();
            numberOfConstraints += nConst*offset.getNumberOfSources();
            numberOfParameters += offset.getNumCalibrationParameters();
        }

        // Define Workspace Size off of DATA
        workSpace.xMax = -std::numeric_limits<double>::max();
        workSpace.yMax = -std::numeric_limits<double>::max();
        workSpace.zMax = -std::numeric_limits<double>::max();
        workSpace.xMin = std::numeric_limits<double>::max();
        workSpace.yMin = std::numeric_limits<double>::max();
        workSpace.zMin = std::numeric_limits<double>::max();
        for( measIT = dataList.begin(); measIT != dataList.end(); measIT ++)
        {
            workSpace.xMax = std::max(workSpace.xMax, measIT->Position.x());
            workSpace.xMin = std::min(workSpace.xMin, measIT->Position.x());

            workSpace.yMax = std::max(workSpace.yMax, measIT->Position.y());
            workSpace.yMin = std::min(workSpace.yMin, measIT->Position.y());

            workSpace.zMax = std::max(workSpace.zMax, measIT->Position.z());
            workSpace.zMin = std::min(workSpace.zMin, measIT->Position.z());
        }

        // Calculate Workspace Info
        pCenter = Eigen::Vector3d(workSpace.xMax+workSpace.xMin,workSpace.yMax+workSpace.yMin,workSpace.zMax+workSpace.zMin)/2.0;

        if( minimumSourceToCenterDistance < 0)
        {
            rMinSq = 0;
            for( measIT = dataList.begin(); measIT != dataList.end(); measIT ++)
            {
                rMinSq = std::max(rMinSq, 1.001*(measIT->Position-pCenter).squaredNorm()); // limit it to no closer than 5% bigger than the wkspace size
            }
        } else
        {
            rMinSq = 1.001 * std::pow(minimumSourceToCenterDistance,2);
        }


        if( maximumSourceToCenterDistance < 0  || (constraint != HEADING_AND_POSITION ))
        {
            rMaxSq = rMinSq;

            if( maximumSourceToCenterDistance > 0 )
                rMaxSq = std::max(std::pow(maximumSourceToCenterDistance,2),rMaxSq);


            int srcNum;
            for(coilIT = coilList.begin(); coilIT != coilList.end(); coilIT ++ )
            {
                ScalorPotential::srcStruct src;
                for( srcNum = 0; srcNum < coilIT->getNumberOfSources(); srcNum ++ )
                {
                    src = coilIT->getSourceStruct(srcNum);
                    rMaxSq = std::max(rMaxSq, 100*(src.srcPosition-pCenter).squaredNorm()); // limit it to 30x the farthest initial guess

                    if( rMaxSq == 0 )
                        rMaxSq = 1;
                }
            }

        } else
        {
            rMaxSq = 0.999 * std::pow(maximumSourceToCenterDistance,2);
        }

        assert(('Electromagnet_Calibration::calibrate: The minimum source to center distance is greator than the maximum allowable source to center distance.', rMaxSq >= rMinSq || constraint == UNIT_HEADING_ONLY ));

        // initialize coefficients with linear least squares
        linearLeastSquareCoeffFit(dataList);


        // preallocate vectors and matricies for solution
        Eigen::MatrixXd J;
        Eigen::MatrixXd JtJ;
        Eigen::VectorXd E;

        // ************ OPTIMIZE WITH NONLINEAR LEAST SQUARES ***************** //


        JtJ.setZero(0,0);
        J.setZero(3*numberOfMeasurements+numberOfConstraints,numberOfParameters);
        E.setZero(0);




        // initialize PHI leave all Lambda's Initialized to Zero
        Eigen::VectorXd states(numberOfParameters), states_last(numberOfParameters), delta_States(numberOfParameters), delta_States_last(numberOfParameters), error_this(numberOfMeasurements*3+numberOfConstraints), delta_error(numberOfMeasurements*3+numberOfConstraints);
        Eigen::VectorXd deltaErrExp(numberOfMeasurements*3+numberOfConstraints);
        Eigen::VectorXd ds_tmp;
        Eigen::VectorXd error_last = error_this;
        Eigen::VectorXd delta_error_last = delta_error;
        Eigen::VectorXd error_last_last = error_last;

        obtainPHI(states);
        applyPHI(states);
        obtainPHI(states);
        states_last = states;

        delta_States.setZero(numberOfParameters);
        delta_States_last.setZero(numberOfParameters);

        packError(error_this, dataList);
        delta_error.setZero(error_this.rows(),error_this.cols());
        error_last = error_this;

        rmsError_this = std::sqrt(error_this.squaredNorm()/error_this.rows());
        rmsError_last = rmsError_this;


        double Lambda = 1;
        double LambdaChangeRange = 1 + sqrtEpsilon;
        double lambdaRangeCount = 0;

        double deltaRmsError = rmsError_this = rmsError_last;

        double predictionFactor = 1.0;
        double percentStateChange = 0.0;
        double percentErrorChange = 0.0;
        int errorIncreaseCounter = 0;



        if( printProgress )
        {
            std::cout << std::setprecision(3) << "Iteration: " << iteration << "   RMS Error: " << rmsError_this << std::endl;// << "   Delta RMS Error: " <<  deltaRmsError  << "\tError Change: " << percentErrorChange << "%" << "\tState Change: " << percentStateChange << "%";
        }

        do {
            converged = true;

            // get updated jacobian
            packErrorJacobian(J, dataList);

            // Remove Zero columns and rescale
            std::vector<int> HtoJ_index;
            std::vector<double> H_col_scale;
            Eigen::MatrixXd H(J.rows(),J.cols());
            unsigned int count = 0;
            for( unsigned int i=0; i<J.cols(); i++ )
            {
                double J_col_norm = J.col(i).norm();
                //cout << "Col: " << i << " Norm: " << J_col_norm << (J_col_norm > std::numeric_limits<double>::epsilon()*10.0?"\tKeep":"\tLoose")<< endl;
                if( J_col_norm > std::numeric_limits<double>::epsilon()*10.0 /*!=0*/ )
                {
                    H.col(count) = J.col(i)/J_col_norm;
                    H_col_scale.push_back(J_col_norm);
                    HtoJ_index.push_back(i);
                    count ++;
                }
            }


            delta_States_last = delta_States;

            int H_rows = H.rows();
            int H_cols = HtoJ_index.size();

            assert( H_rows > 0 && H_cols > 0 );

            JtJ = H.block(0,0,H_rows,H_cols).transpose()*H.block(0,0,H_rows,H_cols);
            JtJ += Lambda*Eigen::MatrixXd::Identity(H_cols,H_cols);

            ds_tmp = -JtJ.fullPivLu().solve(H.block(0,0,H_rows,H_cols).transpose()*error_this);

            //ds_tmp = -JtJ.ldlt().solve(H.block(0,0,H_rows,H_cols).transpose()*error_this);

            delta_States.setZero(delta_States.rows());
            for( unsigned int i=0; i<H_cols; i++ )
            {
                delta_States(HtoJ_index[i]) = ds_tmp(i)/H_col_scale[i];
            }


            states += delta_States;

            delta_error_last = delta_error;
            error_last_last = error_last;
            error_last = error_this;


            applyPHI(states);
            packError(error_this, dataList);


            double rmsError_last_last = rmsError_last;
            rmsError_last = rmsError_this;
            rmsError_this = std::sqrt(error_this.squaredNorm()/error_this.rows());

            delta_error = error_this-error_last;

            deltaRmsError = rmsError_this-rmsError_last;

            deltaErrExp =  J*delta_States + error_last;

            double delta_error_prediction = std::sqrt(deltaErrExp.squaredNorm()/error_this.rows())-rmsError_last;
            predictionFactor = deltaRmsError/( delta_error_prediction );

            delta_States = states-states_last; // just incase BC/Constraint inforcement changed something

            percentErrorChange = deltaRmsError / rmsError_this * 100.0;
            if( std::isnan(percentErrorChange) )
                percentErrorChange = deltaRmsError;

            double normStates = states.norm();
            percentStateChange = delta_States.norm()/normStates * 100.0;
            if( std::isnan(percentStateChange) )
                percentStateChange = sqrt(delta_States.squaredNorm()/delta_States.rows());


            if( deltaRmsError <= 0  )
            {
                errorIncreaseCounter = 1;



                // UPDATE Lambda based on convergance criteria

                // handel boundry conditions that occure at initialization
                if( std::isinf(predictionFactor) || std::isnan(predictionFactor)  )
                {
                    predictionFactor = 1.0-3.0/4.0*LambdaChangeRange;
                }

                if( predictionFactor > 1+LambdaChangeRange/4.0 )
                {
                    double e_this_sqNorm = error_this.squaredNorm();
                    double deltaError = e_this_sqNorm-error_last.squaredNorm();

                    std::vector<double> e,ac,bc,cc;
                    e.push_back(e_this_sqNorm-deltaError);
                    e.push_back(e_this_sqNorm);
                    ac.push_back(0);
                    ac.push_back(1);
                    bc.push_back(0);
                    bc.push_back(1);
                    cc.push_back(1);
                    cc.push_back(1);

                    double factor = 1.0;
                    while (  deltaError < 0 )
                    {
                        factor *= 3.0;
                        states = states_last + factor*delta_States;

                        // update error
                        error_last = error_this;
                        applyPHI(states);
                        packError(error_this, dataList);

                        double e_this_sqNorm2 = error_this.squaredNorm();

                        e.push_back(e_this_sqNorm2);
                        ac.push_back(factor*factor);
                        bc.push_back(factor);
                        cc.push_back(1.0);

                        deltaError = e_this_sqNorm2-e_this_sqNorm;
                        e_this_sqNorm = e_this_sqNorm2;
                    }

                    Eigen::MatrixXd A(e.size(),3);
                    Eigen::VectorXd E(e.size(),1);
                    for( int i=0; i< e.size(); i++ )
                    {
                        A.block(i,0,1,3) << ac[i], bc[i], cc[i];
                        E(i) = e[i];
                    }
                    Eigen::Vector3d coeff = A.householderQr().solve(E);

                    factor = -coeff(1)/(2.0*coeff(0));

                    delta_States *= factor;
                    states = states_last + delta_States;

                    // update error
                    error_last = error_this;
                    applyPHI(states);
                    packError(error_this, dataList);


                    delta_States = states-states_last; // just incase BC/Constraint inforcement changed something

                    error_last = error_last_last;
                    rmsError_this = std::sqrt(error_this.squaredNorm()/error_this.rows());//solutionRMSError();

                    delta_error = error_this-error_last;

                    deltaRmsError = rmsError_this-rmsError_last;

                    lambdaRangeCount ++;

                    percentErrorChange = deltaRmsError / rmsError_this * 100.0;
                    if( std::isnan(percentErrorChange) )
                        percentErrorChange = deltaRmsError;

                    double normStates = states.norm();
                    percentStateChange = delta_States.norm()/normStates * 100.0;
                    if( std::isnan(percentStateChange) )
                        percentStateChange = sqrt(delta_States.squaredNorm()/delta_States.rows());

                }
                else if( std::abs(1.0-predictionFactor) > LambdaChangeRange )
                {
                    // Error change was too small compared to a linear change. Getting close to an over-step.  Be more gradient decent.
                    double factor =  2;
                    Lambda = std::max( Lambda*factor, std::numeric_limits<double>::epsilon() );

                }
                else if( std::abs(1.0-predictionFactor) < LambdaChangeRange/2.0 )
                {
                    // Error change was too close to linear.  Get more aggresive!
                    double factor  =  1/2.0;
                    Lambda = std::max( Lambda*factor, std::numeric_limits<double>::epsilon() );

                    lambdaRangeCount ++;

                }

                if( lambdaRangeCount > 3 )
                {
                    LambdaChangeRange = std::min(1.0, LambdaChangeRange*2);
                    lambdaRangeCount = 0;
                }



                states_last = states;
            }
            else
            {
                errorIncreaseCounter ++;
                // change made things worse!


                Lambda *= 2.0*errorIncreaseCounter;
                //Lambda *= std::abs(deltaRmsError-delta_error_prediction)/Lambda+2.0;//2.0;
                LambdaChangeRange = std::max(LambdaChangeRange*0.75, sqrtEpsilon);


                states = states_last;
                delta_States = delta_States_last;

                error_this = error_last;
                error_last = error_last_last;
                delta_error = delta_error_last;

                rmsError_this = rmsError_last;
                rmsError_last = rmsError_last_last;

                deltaRmsError = rmsError_this-rmsError_last;

                delta_States_last.setZero(delta_States_last.rows(),delta_States_last.cols());
                lambdaRangeCount = 0;

                converged = false;
            }

            iteration ++;

            assert( !std::isinf(Lambda) && !std::isnan(Lambda));
            assert( !std::isinf(rmsError_this) && !std::isnan(rmsError_this));



            converged = converged && (
                        iteration >= maxAttempts ||
                        rmsError_this <converganceTolerance ||
                        (std::abs(deltaRmsError)*(1+Lambda) < converganceTolerance && delta_error.norm()!=0) ||
                        (std::abs(percentErrorChange)*(1+Lambda)/100.0 < converganceTolerance && percentErrorChange >= std::numeric_limits<double>::epsilon()*rmsError_this ) ||
                        (std::abs(percentStateChange)*(1+Lambda)/100.0 < converganceTolerance && percentStateChange >= std::numeric_limits<double>::epsilon()*normStates ) );


            // make sure we have multiple converged runs before declaring done
            if( deltaRmsError < 0 )
            {
                convergedIterations = converged*(convergedIterations+1);
                converged = converged && (convergedIterations >= numberOfConvergedIterations);
            }

            if( iteration == 1 && maxAttempts > 1 )
                converged = false;

            if( Lambda > 10000 )
                converged = true;


            if( printProgress && iteration%100 == 1)
            {
                std::cout <</* "Conv: " << (converged?"YES   ":"NO   ") << */std::setprecision(3) << "Iteration: " << iteration << "\tRMS Error: " << 1000*std::sqrt(error_this.head(3*numberOfMeasurements).squaredNorm()/(3*numberOfMeasurements)) << " mT\tDelta RMS Error: " <<  deltaRmsError*1000  << " mT\t  Error Change: " << percentErrorChange << "%" << "\tState Change: " << percentStateChange << "%";
                std::cout << std::endl;
            }


        }while( !converged );

        if( printProgress )
        {
            std::cout << std::setprecision(3) << "Iteration: " << iteration << "\tRMS Error: " << 1000*std::sqrt(error_this.head(3*numberOfMeasurements).squaredNorm()/(3*numberOfMeasurements)) << " mT\tDelta RMS Error: " <<  deltaRmsError*1000  << " mT\tPercent Error Change: " << percentErrorChange << "%" << "\tPercent State Change: " << percentStateChange << "%\tTolerance: " << converganceTolerance << "\tDONE!"<<std::endl<<std::endl;
        }


        applyPHI(states);

        if( constraint != UNIT_HEADING_ONLY && ! checkSourcePositions( printProgress ) )
        {
            if( maximumSourceToCenterDistance > 0 )
                rMaxSq = std::max(rMinSq*1.001, std::pow(maximumSourceToCenterDistance,2));

            minRadIsActive = true;
            nConst = 2;
            posWeight *= 2;
            converged = false;

            if( printProgress )
                cout << "Resolving to move sources out of workspace" << endl;
        }

    }while(!converged);


    // remove any residual errors
    // initialize coefficients with linear least squares
    linearLeastSquareCoeffFit(dataList);


    if( printStats_ )
        printStats( dataList);

    return;
}

void ElectromagnetCalibration::linearLeastSquareCoeffFit(const std::vector<MagneticMeasurement> &dataList)
{
    // ************ INITIALIZE COEFFICIENTS WITH LINEAR LEAST SQUARES ***************** //
    Eigen::VectorXd E(numberOfMeasurements*3);
    Eigen::MatrixXd JtJ(numberOfParameters - numberOfSources*6, numberOfParameters-numberOfSources*6);
    Eigen::MatrixXd J(3*numberOfMeasurements, numberOfParameters - numberOfSources*6);

    // Pack Error
    unsigned int phiRowNum = 0;
    unsigned int coilNum;
    std::vector<ScalorPotential>::iterator coilIT;
    std::vector<MagneticMeasurement>::const_iterator measIT;
    for( coilIT = coilList.begin(), coilNum = 0;
         coilIT != coilList.end();
         coilIT ++, coilNum ++ )
    {
        // for each source
        int srcNum;
        std::vector<ScalorPotential::srcStruct>::iterator srcIT;
        for( srcIT = coilIT->srcList.begin(), srcNum = 0;
             srcIT != coilIT->srcList.end();
             srcIT ++, srcNum++ )
        {
            int numParam = coilIT->getNumCalibrationParameters( srcNum );

            int measColInd;
            for( measIT = dataList.begin(), measColInd = 0;
                 measIT != dataList.end();
                 measIT ++, measColInd+=3)
            {
                Eigen::MatrixXd Ft,Hessian;
                Ft.setZero(numParam,3);
                Hessian.setZero(numParam,numParam);
                coilIT->packJacobians(srcNum, measIT->AppliedCurrentVector(coilNum), measIT->Position, -measIT->Field, Ft, Hessian);


                J.block(measColInd, phiRowNum, 3, numParam-6 ) = Ft.topRows(numParam-6).transpose();

                E.segment(measColInd,3) = -measIT->Field;
            }
            phiRowNum += numParam-6;
        }
    }

    if( useOffset() )
    {
        // for each offset source
        int srcNum;
        std::vector<ScalorPotential::srcStruct>::iterator srcIT;
        for( srcIT = offset.srcList.begin(), srcNum = 0;
             srcIT != offset.srcList.end();
             srcIT ++, srcNum++ )
        {
            int numParam = offset.getNumCalibrationParameters( srcNum );

            int measColInd;
            for( measIT = dataList.begin(), measColInd = 0;
                 measIT != dataList.end();
                 measIT ++,  measColInd+=3)
            {
                Eigen::MatrixXd Ft,Hessian;
                Ft.setZero(numParam,3);
                Hessian.setZero(numParam,numParam);
                offset.packJacobians(srcNum, 1.0, measIT->Position, -measIT->Field, Ft, Hessian);


                J.block(measColInd, phiRowNum, 3, numParam-6) = Ft.topRows(numParam-6).transpose();

                E.segment(measColInd, 3) = (-measIT->Field);
            }
            phiRowNum += numParam-6;
        }
    }

    JtJ = J.transpose()*J;
    Eigen::VectorXd PHI = -JtJ.fullPivLu().solve(J.transpose()*E);

    // apply initial coefficient fit to parameters
    phiRowNum = 0;
    for( coilIT = coilList.begin(), coilNum = 0;
         coilIT != coilList.end();
         coilIT ++, coilNum ++ )
    {
        // for each source
        int srcNum;
        std::vector<ScalorPotential::srcStruct>::iterator srcIT;
        for( srcIT = coilIT->srcList.begin(), srcNum = 0;
             srcIT != coilIT->srcList.end();
             srcIT ++, srcNum++ )
        {

            std::vector<ScalorPotential::srcCoeff>::iterator coeffIT;
            for( coeffIT = srcIT->A_Coeff.begin(); coeffIT != srcIT->A_Coeff.end(); coeffIT ++, phiRowNum++ )
            {
                coeffIT->coeff = PHI(phiRowNum);
            }
            for( coeffIT = srcIT->B_Coeff.begin(); coeffIT != srcIT->B_Coeff.end(); coeffIT ++, phiRowNum++ )
            {
                coeffIT->coeff = PHI(phiRowNum);
            }

        }
    }

    if( useOffset() )
    {
        // for each offset source
        int srcNum;
        std::vector<ScalorPotential::srcStruct>::iterator srcIT;
        for( srcIT = offset.srcList.begin(), srcNum = 0;
             srcIT != offset.srcList.end();
             srcIT ++, srcNum++ )
        {
            std::vector<ScalorPotential::srcCoeff>::iterator coeffIT;
            for( coeffIT = srcIT->A_Coeff.begin(); coeffIT != srcIT->A_Coeff.end(); coeffIT ++, phiRowNum++ )
            {
                coeffIT->coeff = PHI(phiRowNum);
            }
            for( coeffIT = srcIT->B_Coeff.begin(); coeffIT != srcIT->B_Coeff.end(); coeffIT ++, phiRowNum++ )
            {
                coeffIT->coeff = PHI(phiRowNum);
            }
        }
    }

}

void ElectromagnetCalibration::applyPHI(const Eigen::VectorXd& PHI )
{
    // initialize PHI leave all Lambda's Initialized to Zero
    unsigned int paramCount = 0;
    for( std::vector<ScalorPotential>::iterator coilIT = coilList.begin(); coilIT != coilList.end(); coilIT ++)
    {
        int numParam = coilIT->getNumCalibrationParameters();
        coilIT->unpackCalibrationState( PHI.segment(paramCount,numParam) );
        paramCount += numParam;
    }
    if( useOffset() )
    {
        int numParam = offset.getNumCalibrationParameters();
        offset.unpackCalibrationState( PHI.segment(paramCount,numParam) );
        paramCount += numParam;
    }
}

void ElectromagnetCalibration::obtainPHI( Eigen::VectorXd& PHI )
{
    // initialize PHI leave all Lambda's Initialized to Zero
    unsigned int paramCount = 0;
    for( std::vector<ScalorPotential>::iterator coilIT = coilList.begin(); coilIT != coilList.end(); coilIT ++)
    {
        int numParam = coilIT->getNumCalibrationParameters();
        coilIT->packCalibrationState( PHI.segment(paramCount,numParam) );
        paramCount += numParam;
    }
    if( useOffset() )
    {
        int numParam = offset.getNumCalibrationParameters();
        offset.packCalibrationState( PHI.segment(paramCount,numParam) );
        paramCount += numParam;
    }
}

void ElectromagnetCalibration::packError( Eigen::VectorXd& error, const std::vector<MagneticMeasurement> & dataList )
{
    // update measurement error
    assert(('Error Vector is the Wrong Size', error.rows() == 3*numberOfMeasurements + numberOfConstraints));

    int measInd;
    std::vector<MagneticMeasurement>::const_iterator measIT;
    for( measIT = dataList.begin(), measInd = 0;
         measIT != dataList.end();
         measIT ++, measInd += 3 )
    {
        assert( ('ElectromagnetCalibration::packError:  Measurement Index Out of Bounds', measInd+2 < 3*dataList.size()) );
        error.segment(measInd,3) = fieldAtPoint(measIT->AppliedCurrentVector,measIT->Position)-measIT->Field;
    }

    if( nConst > 0 )
    {
        // update Constraint Error
        std::vector<ScalorPotential>::iterator coilIT;
        for( coilIT = coilList.begin(); coilIT != coilList.end(); coilIT ++ )
        {
            // for each source
            std::vector<ScalorPotential::srcStruct>::iterator srcIT;
            for( srcIT = coilIT->srcList.begin();
                 srcIT != coilIT->srcList.end();
                 srcIT ++, measInd+=nConst)
            {

                assert( ('ElectromagnetCalibration::packError:  Measurement Index Out of Bounds', measInd < error.rows()) );
                error(measInd) = srcIT->srcDirection.squaredNorm()-1;

                if( nConst == 2 )
                {
                    // srcPosition Constraint Addition
                    double rSqSrc = (srcIT->srcPosition - pCenter).squaredNorm();
                    double rSqRatio = 1;
                    double rSqW = 1;

                    if( minRadIsActive && rSqSrc <= rMinSq && srcIT->A_Coeff.size() == 0 )
                    {
                        rSqRatio = rSqSrc/rMinSq;
                        rSqW = 1.001*rMinSq;
                    }
                    else if(rSqSrc >= rMaxSq )
                    {
                        rSqRatio = rSqSrc/rMaxSq;
                        rSqW = 0.999*rMaxSq;
                    }

                    error(measInd+1) = posWeight*(rSqRatio - 1);

                }
            }
        }

        if( useOffset() )
        {
            // for each offset source
            std::vector<ScalorPotential::srcStruct>::iterator srcIT;
            for( srcIT = offset.srcList.begin(); srcIT != offset.srcList.end(); srcIT ++, measInd += nConst)
            {

                assert( ('ElectromagnetCalibration::packError:  Measurement Index Out of Bounds', measInd < error.rows()) );
                error(measInd) = srcIT->srcDirection.squaredNorm()-1;


                if( nConst == 2)
                {
                    // srcPosition Constraint Addition
                    double rSqSrc = (srcIT->srcPosition - pCenter).squaredNorm();
                    double rSqRatio = 1;
                    double rSqW = 1;

                    if( minRadIsActive && rSqSrc <= rMinSq && srcIT->A_Coeff.size() == 0 )
                    {
                        rSqRatio = rSqSrc/rMinSq;
                        rSqW = 1.001*rMinSq;
                    }
                    else if(rSqSrc >= rMaxSq )
                    {
                        rSqRatio = rSqSrc/rMaxSq;
                        rSqW = 0.999*rMaxSq;
                    }

                    error(measInd+1) = posWeight*(rSqRatio - 1);


                }
            }

        }
    }
}

void ElectromagnetCalibration::packErrorJacobian( Eigen::MatrixXd& J, const std::vector<MagneticMeasurement> & dataList  )
{  
    // Pack Jacobian
    unsigned int phiRowNum = 0;
    unsigned int coilNum;
    std::vector<ScalorPotential>::iterator coilIT;
    for( coilIT = coilList.begin(), coilNum = 0;
         coilIT != coilList.end();
         coilIT ++, coilNum ++ )
    {
        // for each source
        int srcNum;
        std::vector<ScalorPotential::srcStruct>::iterator srcIT;
        for( srcIT = coilIT->srcList.begin(), srcNum = 0;
             srcIT != coilIT->srcList.end();
             srcIT ++, srcNum++ )
        {
            int numParam = coilIT->getNumCalibrationParameters( srcNum );

            int measColInd;
            std::vector<MagneticMeasurement>::const_iterator measIT;
            for(measIT = dataList.begin(), measColInd = 0;
                measIT != dataList.end();
                measIT ++, measColInd+=3)
            {
                Eigen::MatrixXd Ft,Hessian;
                Ft.setZero(numParam,3);
                Hessian.setZero(numParam,numParam);
                coilIT->packJacobians(srcNum, measIT->AppliedCurrentVector(coilNum), measIT->Position, Eigen::Vector3d(0,0,0), Ft, Hessian);


                J.block(measColInd, phiRowNum, 3, numParam ) = Ft.transpose();
            }
            phiRowNum += numParam;
        }
    }

    if( useOffset() )
    {
        // for each offset source
        int srcNum;
        std::vector<ScalorPotential::srcStruct>::iterator srcIT;
        for( srcIT = offset.srcList.begin(), srcNum = 0;
             srcIT != offset.srcList.end();
             srcIT ++, srcNum++ )
        {
            int numParam = offset.getNumCalibrationParameters( srcNum );

            int measColInd;
            std::vector<MagneticMeasurement>::const_iterator measIT;
            for(measIT = dataList.begin(), measColInd = 0;
                measIT != dataList.end();
                measIT ++, measColInd+=3)
            {
                Eigen::MatrixXd Ft,Hessian;
                Ft.setZero(numParam,3);
                Hessian.setZero(numParam,numParam);
                offset.packJacobians(srcNum, 1.0, measIT->Position, Eigen::Vector3d(0,0,0), Ft, Hessian);

                J.block(measColInd, phiRowNum, 3, numParam ) = Ft.transpose();;
            }
            phiRowNum += numParam;
        }
    }


    if( nConst > 0 )
    {
        // Add In Constraint Terms
        int constraintNum = 0;
        phiRowNum = 0;
        int measColInd = numberOfMeasurements*3;
        for( coilIT = coilList.begin(); coilIT != coilList.end(); coilIT ++ )
        {
            // for each source
            int srcNum;
            std::vector<ScalorPotential::srcStruct>::iterator srcIT;
            for( srcIT = coilIT->srcList.begin(), srcNum = 0; srcIT != coilIT->srcList.end();
                 srcIT ++, srcNum++, measColInd+=nConst)
            {
                int numParam = coilIT->getNumCalibrationParameters( srcNum );
                J.block(measColInd, phiRowNum+numParam-6, 1, 3) = 2*srcIT->srcDirection.transpose();

                if( nConst ==2 )
                {
                    // srcPosition Constraint Addition
                    double rSqSrc = (srcIT->srcPosition - pCenter).squaredNorm();
                    double rSqRatio = 1;
                    double rSqW = 1;

                    if( minRadIsActive && rSqSrc <= rMinSq && srcIT->A_Coeff.size() == 0 )
                    {
                        rSqRatio = rSqSrc/rMinSq;
                        rSqW = 1.001*rMinSq;
                    }
                    else if(rSqSrc >= rMaxSq )
                    {
                        rSqRatio = rSqSrc/rMaxSq;
                        rSqW = 0.999*rMaxSq;
                    }

                    J.block(measColInd+1, phiRowNum+numParam-3, 1, 3) = 2.0*posWeight*(srcIT->srcPosition - pCenter).transpose()/rSqW;

                }

                phiRowNum += numParam;
            }
        }

        if( useOffset() )
        {
            // for each offset source
            int srcNum;
            std::vector<ScalorPotential::srcStruct>::iterator srcIT;
            for( srcIT = offset.srcList.begin(), srcNum = 0; srcIT != offset.srcList.end();
                 srcIT ++, srcNum++, measColInd += nConst)
            {
                int numParam = offset.getNumCalibrationParameters( srcNum );

                J.block(measColInd, phiRowNum+numParam-6, 1, 3) = 2*srcIT->srcDirection.transpose();

                if( nConst == 2)
                {
                    // srcPosition Constraint Addition
                    double rSqSrc = (srcIT->srcPosition - pCenter).squaredNorm();
                    double rSqRatio = 1;
                    double rSqW = 1;

                    if( minRadIsActive && rSqSrc <= rMinSq && srcIT->A_Coeff.size() == 0 )
                    {
                        rSqRatio = rSqSrc/rMinSq;
                        rSqW = 1.001*rMinSq;
                    }
                    else if(rSqSrc >= rMaxSq )
                    {
                        rSqRatio = rSqSrc/rMaxSq;
                        rSqW = 0.999*rMaxSq;
                    }

                    J.block(measColInd+1, phiRowNum+numParam-3, 1, 3) = 2.0*posWeight*(srcIT->srcPosition - pCenter).transpose()/rSqW;

                }

                phiRowNum += numParam;
            }
        }
    }
}


void ElectromagnetCalibration::printStats( const std::vector< MagneticMeasurement >& dataList) const
{


    double percentError = 0;
    double SS_res = 0;
    double SS_total = 0;
    double SS_field = 0;

    Eigen::Vector3d avgFieldError(0,0,0);
    Eigen::Vector3d fieldTmp;
    Eigen::Vector3d avgField(Eigen::Vector3d::Zero());

    std::vector<MagneticMeasurement>::const_iterator dataIT;

    for( dataIT  = dataList.begin();
         dataIT != dataList.end();
         dataIT++ )
    {
        SS_field += dataIT->Field.squaredNorm();

    }
    double rmsField = std::sqrt(SS_field/(double)dataList.size()/3.0);

    double averageFielddMag = 0;
    for( dataIT  = dataList.begin();
         dataIT != dataList.end();
         dataIT++ )
    {
        avgField += dataIT->Field;
        averageFielddMag += dataIT->Field.norm();
        fieldTmp = fieldAtPoint(dataIT->AppliedCurrentVector, dataIT->Position );
        avgFieldError += (dataIT->Field - fieldTmp);
        percentError += sqrt((dataIT->Field - fieldTmp).squaredNorm()/3.0)/rmsField;//dataIT->Field.norm();
        SS_res += (dataIT->Field - fieldTmp).squaredNorm();
        SS_total += dataIT->Field.squaredNorm();
    }

    avgField /= (double) dataList.size();
    averageFielddMag /= (double) dataList.size();
    avgFieldError /= (double) dataList.size();
    Eigen::Matrix3d errorVariance; errorVariance.setZero(3,3);


    Eigen::Matrix3d fieldVar(Eigen::Matrix3d::Zero());
    for( dataIT  = dataList.begin();
         dataIT != dataList.end();
         dataIT++ )
    {


        fieldTmp = fieldAtPoint(dataIT->AppliedCurrentVector, dataIT->Position );
        Eigen::Vector3d fieldError = (dataIT->Field - fieldTmp);
        Eigen::Vector3d errorDiff = fieldError-avgFieldError;
        errorVariance += errorDiff*errorDiff.transpose();

        errorDiff = dataIT->Field-avgField;
        fieldVar += errorDiff*errorDiff.transpose();
    }

    fieldVar /= ((double) dataList.size()-1);
    errorVariance /= ((double) dataList.size()-1);

    percentError /= (double) dataList.size()/100.0;
    double rmsError_out = std::sqrt(SS_res/(double)dataList.size()/3.0);
    double rmsError_normalized = rmsError_out/std::sqrt(SS_field/(double)dataList.size()/3.0);
    double R_squared_out = 1- SS_res/SS_total;



    Eigen::JacobiSVD<Eigen::Matrix3d> varSVD(errorVariance,Eigen::ComputeFullU|Eigen::ComputeFullV);
    Eigen::Matrix3d S = varSVD.singularValues().asDiagonal();
    for( int i=0;i<3;i++ )
        S(i,i) = sqrt(S(i,i));
    Eigen::Matrix3d stdFieldErr = varSVD.matrixU()*S*varSVD.matrixV().transpose();

    cout << getName() << endl
         << "  Percent Error:\t" << setprecision(4) << percentError << "%" << endl
         << "  RMS Error:\t" << setprecision(4) << rmsError_out << endl
         << "  RMS Error/RMS Field: " << setprecision(4) << rmsError_normalized*100 <<"%"<< endl
         << "  R^2:\t\t" << setprecision(9) << R_squared_out << endl
         << "  Average Field Magnitude: " << setprecision(4) << averageFielddMag*1000 << " mT" << endl
         << "  Average Field Error:\t(" << setprecision(4) << (avgFieldError.transpose()*1000) << ") mT" << endl
         << "  Average Field Error:\t(" << setprecision(4) << (avgFieldError.transpose()/averageFielddMag)*100 << ") %" << endl
         << "  \nField Error Variance:\n" << "---------------\n" << setprecision(4) << avgFieldError << endl << "---------------" << endl
         << "  \nField Error stdev [mT]:\n" << "---------------\n" << setprecision(4) << stdFieldErr*1000 << endl << "---------------" << endl
         << "  \nField Error stdev [%]:\n" << "---------------\n" << setprecision(4) << stdFieldErr/averageFielddMag*100 << endl << "---------------" << endl;

    cout << endl;

    return;
}
