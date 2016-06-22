#include "scalorPotential.h"
#include <iomanip>

ScalorPotentialState::ScalorPotentialState()
{
    value = 0;
    firstSpatialDerivative.setZero(3,1);
    secondSpatialDerivative.setZero(3,3);
    thirdSpatialDerivative.setZero(5,3);

    firstSpatialDerivative_SourceHeadingDerivative.clear();
    secondSpatialDerivative_SourceHeadingDerivative.clear();

    firstSpatialDerivative_SourcePositionDerivative.clear();
    secondSpatialDerivative_SourcePositionDerivative.clear();

}

ScalorPotentialCalibrationJacobians::ScalorPotentialCalibrationJacobians()
{
    firstSpatialDerivative.setZero(3,1); ///< Field, the gradient of the potential.
    firstSpatialDerivative_A_CoeffDerivative.setZero(0,0); ///< How the field changes with the A coefficients
    firstSpatialDerivative_B_CoeffDerivative.setZero(0,0); ///< How the field changes with the B coefficients

    firstSpatialDerivative_SourcePositionDerivative.setZero(3,3); ///< Field spatial gradient
    firstSpatialDerivative_SourceHeadingDerivative.setZero(3,3); ///< How the field changes with the source heading
    firstSpatialDerivative_secondSourceHeadingDerivative.setZero(9,3); ///< the second deritive of field with heading (d(BX)dz; d(BY/dz; d(BZ)/dz)


    firstSpatialDerivative_dA_dHeading.setZero(0,0);
    firstSpatialDerivative_dB_dHeading.setZero(0,0);
    firstSpatialDerivative_dA_dPosition.setZero(0,0);
    firstSpatialDerivative_dB_dPosition.setZero(0,0);
}

ScalorPotential::srcCoeff::srcCoeff():
    order(0),
    coeff(0)
{
    return;
}

ScalorPotential::srcCoeff::srcCoeff(double value, unsigned int order):
    order(order),
    coeff(value)
{
    return;
}

ScalorPotential::srcStruct::srcStruct():
    A_Coeff(std::vector<srcCoeff>(0)),
    B_Coeff(std::vector<srcCoeff>(0)),
    srcPosition(0,0,0),
    srcDirection(0,0,1)
{
    return;
}


unsigned int ScalorPotential::srcStruct::getMaxOrder_A_Coeff() const
{
    unsigned int maxCoeffOrder = 0;
    for( unsigned int i=0; i<A_Coeff.size(); i++ )
        maxCoeffOrder = std::max(maxCoeffOrder,A_Coeff[i].order );

    return maxCoeffOrder;
}

unsigned int ScalorPotential::srcStruct::getMaxOrder_B_Coeff() const
{
    unsigned int maxCoeffOrder = 0;
    for( unsigned int i=0; i<B_Coeff.size(); i++ )
        maxCoeffOrder = std::max(maxCoeffOrder,B_Coeff[i].order );

    return maxCoeffOrder;
}


ScalorPotential::ScalorPotential():
    srcList(std::vector<srcStruct>(0))
{
    return;
}

///
/// \brief the constructor from a coil list.
/// \param coilList The list of coils and their respective sources.
/// \param dc_field_offset The dc offset field, if any.
///
ScalorPotential::ScalorPotential(const std::vector<srcStruct>& srcList_ ):
    srcList(srcList_)
{

    int numA = 0;
    int numB = 0;
    int numP = srcList.size()*3;
    int numZ = srcList.size()*3;

    for( unsigned int i=0; i<srcList.size(); i++ )
    {
        numA += srcList[i].A_Coeff.size();
        numB += srcList[i].B_Coeff.size();
    }

    numCalParameters =  numA + numB + numP + numZ;
}

// This function returns the field, gradient, gradientJacobian, and field/gradient current jacobain.  It is more efficient than requesting them seporately
ScalorPotentialState ScalorPotential::getState(const Eigen::Vector3d& position, int sourceNumber ) const
{

    ScalorPotentialState returnState;

    if( sourceNumber > 0 && sourceNumber < (int)srcList.size() )
    {
        srcFieldGradient(position,srcList[sourceNumber],returnState);

    }else if( sourceNumber == -1 )
    {
        for( unsigned int src = 0; src < srcList.size(); src++ )
        {
            srcFieldGradient(position, srcList[src], returnState );
        }
    }

    return returnState;

}


// These two functions convert a vector packing of the gradient to a matrix packing and vice versa.
Eigen::Matrix3d ScalorPotential::remapSecondDerivativeVec(const Vector5d& gradVector )
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

Vector5d ScalorPotential::remapSecondDerivativeMat( const Eigen::Matrix3d& gradMatrix)
{
    Vector5d gradVec;
    gradVec.topRows<3>() = gradMatrix.leftCols<1>();
    gradVec(3) = gradMatrix(1,1);
    gradVec(4) = gradMatrix(2,1);
    return gradVec;
}

///@brief returns the number of sources for the given coil
unsigned int ScalorPotential::getNumberOfSources( ) const
{
    return srcList.size();
}

ScalorPotential::srcStruct ScalorPotential::getSourceStruct(unsigned int sourceNumber) const
{
    if( sourceNumber < getNumberOfSources() )
        return srcList[sourceNumber];
    else
        return srcStruct();
}

void ScalorPotential::setSourceStruct(unsigned int sourceNumber, const srcStruct& newSrc)
{

    assert(('ScalorPotential::setSourceStruct: newSrc must have a direction of length greator than zero', newSrc.srcDirection.norm() != 0 ));

    if( sourceNumber < getNumberOfSources() )
        srcList[sourceNumber] = newSrc;
    else
    {
        while( sourceNumber !=0 && srcList.size()< sourceNumber-1)
        {
            srcList.push_back(srcStruct());
        }
        srcList.push_back(newSrc);
        srcList.back().srcDirection.normalize();
    }

    int numA = 0;
    int numB = 0;
    int numP = srcList.size()*3;
    int numZ = srcList.size()*3;

    for( int i=0; i<srcList.size(); i++ )
    {
        numA += srcList[i].A_Coeff.size();
        numB += srcList[i].B_Coeff.size();
    }

    numCalParameters =  numA + numB + numP + numZ;
}

void ScalorPotential::removeSourceStruct(unsigned int sourceNumber)
{
    if( sourceNumber < getNumberOfSources() )
    {
        std::vector<srcStruct>::iterator it = srcList.begin();
        for( unsigned int i=0;i=sourceNumber; i++ ) it++;

        srcList.erase(it);
    }

    int numA = 0;
    int numB = 0;
    int numP = srcList.size()*3;
    int numZ = srcList.size()*3;

    for( int i=0; i<srcList.size(); i++ )
    {
        numA += srcList[i].A_Coeff.size();
        numB += srcList[i].B_Coeff.size();
    }

    numCalParameters =  numA + numB + numP + numZ;
}


double ScalorPotential::getValue( const Eigen::Vector3d& position, int sourceNumber ) const
{

    // field Coefficients
    double value = 0;


    int srcStart = sourceNumber;
    int srcEnd = sourceNumber+1;
    if(sourceNumber == -1 )
    {
        srcStart = 0;
        srcEnd = srcList.size();
    } else if( sourceNumber >= srcList.size() )
        return 0;

    for( int srcNum=srcStart; srcNum<srcEnd; srcNum++ )
    {
        Eigen::Vector3d p( position - srcList[srcNum].srcPosition );
        double pMag = p.norm();

        if( pMag > 0 )
        {
            p /= pMag;
        }

        Eigen::Vector3d z( srcList[srcNum].srcDirection );
        if( z.norm() == 0 )
            z(3) = 1;

        z.normalize();

        // define Cosine Theta
        double cTh = p.dot(z);

        std::vector<srcCoeff>::const_iterator coeffIT;

        for( coeffIT = srcList[srcNum].A_Coeff.begin();
             coeffIT !=srcList[srcNum].A_Coeff.end();
             coeffIT ++ )
        {
            double n = coeffIT->order;
            double Pn   = LegandrePolynomial(cTh, n, 0); // Legandre Polynomial at cos(theta)

            value += coeffIT->coeff * std::pow(pMag,n)*Pn;

        }

        if( pMag > 0 )
        {
            for( coeffIT = srcList[srcNum].B_Coeff.begin();
                 coeffIT !=srcList[srcNum].B_Coeff.end();
                 coeffIT ++ )
            {
                double n = coeffIT->order;
                double Pn   = LegandrePolynomial(cTh, n, 0); // Legandre Polynomial at cos(theta)
                value += coeffIT->coeff * 1.0/std::pow(pMag,n+1.0)*Pn;
            }


        }
    }

    return value;

}



Eigen::Vector3d ScalorPotential::getGradient( const Eigen::Vector3d& position, int sourceNumber ) const
{

    Eigen::Vector3d field(0,0,0);


    int srcStart = sourceNumber;
    int srcEnd = sourceNumber+1;
    if(sourceNumber == -1 )
    {
        srcStart = 0;
        srcEnd = srcList.size();
    }else if( sourceNumber >= srcList.size() )
        return Eigen::Vector3d::Zero();

    for( int srcNum=srcStart; srcNum<srcEnd; srcNum++ )
    {
        Eigen::Vector3d p( position - srcList[srcNum].srcPosition );
        double pMag = p.norm();

        if( pMag > 0 )
        {
            p /= pMag;
        }

        Eigen::Vector3d z( srcList[srcNum].srcDirection );
        if( z.norm() == 0 )
            z(3) = 1;

        z.normalize();

        // define Cosine Theta
        double cTh = p.dot(z);

        // field Coefficients
        double field_r_mult=0;
        double field_z_mult=0;

        std::vector<srcCoeff>::const_iterator coeffIT;

        for( coeffIT = srcList[srcNum].A_Coeff.begin();
             coeffIT !=srcList[srcNum].A_Coeff.end();
             coeffIT ++ )
        {
            double n = coeffIT->order;
            double Pn   = LegandrePolynomial(cTh, n, 0); // Legandre Polynomial at cos(theta)
            double Pn_1 = LegandrePolynomial(cTh, n, 1); // First Derivative of Legandre Polynomial at cos(theta)

            if( n==0 )
                continue;


            double c1 = 0;
            double c2 = 0;

            double tmp =  coeffIT->coeff;

            if( n > 1 )
                tmp *= pow(pMag,n-1);

            c2 += tmp;

            tmp *= n;
            c1 += tmp;

            // field terms
            field_r_mult += c2*cTh*Pn_1-c1*Pn;
            field_z_mult -= c2*Pn_1;
        }

        if( pMag > 0 )
        {
            for( coeffIT = srcList[srcNum].B_Coeff.begin();
                 coeffIT !=srcList[srcNum].B_Coeff.end();
                 coeffIT ++ )
            {
                double n = coeffIT->order;
                double Pn   = LegandrePolynomial(cTh, n, 0); // Legandre Polynomial at cos(theta)
                double Pn_1 = LegandrePolynomial(cTh, n, 1); // First Derivative of Legandre Polynomial at cos(theta)


                double c1 = 0;
                double c2 = 0;

                double tmp =  coeffIT->coeff/pow(pMag,n+2);
                c2 += tmp;

                tmp *= (n+1);
                c1 += -1.0*tmp;

                // field terms
                field_r_mult += c2*cTh*Pn_1-c1*Pn;
                field_z_mult -= c2*Pn_1;

            }
        }

        field+= field_r_mult*p + field_z_mult*z;
    }



    // calculate field
    return field;
}


// This fuction returns the field and gradient (stacked) assuming a current of 1 at the position specified in the work space
void ScalorPotential::srcFieldGradient(const Eigen::Vector3d& position, const srcStruct& src, ScalorPotentialState& currentState )
{

    Eigen::Vector3d p( position - src.srcPosition );
    double pMag = p.norm();
    double pMagSq = pMag*pMag;

    if( pMag > 0 )
    {
        p /= pMag;
    }

    Eigen::Vector3d z( src.srcDirection );
    if( z.norm() == 0 )
        z(3) = 1;

    z.normalize();


    Eigen::Matrix3d pzt_zpt( p*z.transpose()+z*p.transpose() );
    Eigen::Matrix3d ppt( p*p.transpose() );// = p*p.transpose(); // p*p'
    Eigen::Matrix3d zzt(z *z.transpose() );// = src.srcDirection *src.srcDirection.transpose(); // z*z'
    Eigen::Matrix3d rzt = p*z.transpose();
    Eigen::Matrix3d zrt = z*p.transpose();

    // some constant vectors and matricies
    Eigen::Vector3d X(1,0,0);
    Eigen::Vector3d Y(0,1,0);
    Eigen::Matrix3d I(Eigen::Matrix3d::Identity());


    // define Cosine Theta
    double cTh = p.dot(z);
    double cThsq = cTh*cTh;

    // field Coefficients
    double field_r_mult=0;
    double field_z_mult=0;

    // gradient Coefficients
    double grad_I_mult=0;
    double grad_rr_mult=0;
    double grad_zz_mult=0;
    double grad_rzzr_mult=0;

    // field change with source heading
    double dBdz_grz_rrt_mult = 0;
    double dBdz_zzt_mult = 0;
    double dBdz_zrt_mult = 0;
    double dBdz_I_mult = 0 ;

    // gradient Jacobian Coefficients
    double gradJacob_zmmz_mult = 0;
    double gradJacob_rmmr_mult = 0;
    double gradJacob_zrrz_zm_mult = 0;
    double gradJacob_zrrz_rm_mult = 0;
    double gradJacob_I_zm_mult = 0;
    double gradJacob_I_rm_mult = 0;
    double gradJacob_zz_zm_mult = 0;
    double gradJacob_zz_rm_mult = 0;
    double gradJacob_rr_zm_mult = 0;
    double gradJacob_rr_rm_mult = 0;

    // gradZhatDerJacob
    double gradZJacob_I_rm_mult  = 0;
    double gradZJacob_I_zm_mult  = 0;
    double gradZJacob_rr_rm_mult = 0;
    double gradZJacob_rr_zm_mult = 0;
    double gradZJacob_zz_rm_mult = 0;
    double gradZJacob_zz_zm_mult = 0;
    double gradZJacob_rz_rm_mult = 0;
    double gradZJacob_rz_zm_mult = 0;
    double gradZJacob_zr_zm_mult = 0;
    double gradZJacob_zr_rm_mult = 0;
    double gradZJacob_mrrm_mult = 0;
    double gradZJacob_mz_mult = 0;
    double gradZJacob_zm_mult = 0;



    std::vector<srcCoeff>::const_iterator coeffIT;

    for( coeffIT = src.A_Coeff.begin();
         coeffIT !=src.A_Coeff.end();
         coeffIT ++ )
    {
        double n = coeffIT->order;
        double Pn   = LegandrePolynomial(cTh, n, 0); // Legandre Polynomial at cos(theta)
        double Pn_1 = LegandrePolynomial(cTh, n, 1); // First Derivative of Legandre Polynomial at cos(theta)
        double Pn_2 = LegandrePolynomial(cTh, n, 2); // Second Derivative of Legandre Polynomial at cos(theta)
        double Pn_3 = LegandrePolynomial(cTh, n, 3); // Third Derivative of Legandre Polynomial at cos(theta)


        currentState.value += coeffIT->coeff * std::pow(pMag,n)*Pn;

        if( n==0 )
            continue;


        double c1 = 0;
        double c2 = 0;
        double c3 = 0;
        double c4 = 0;

        double tmp =  coeffIT->coeff;

        if( n > 1 )
            tmp *= pow(pMag,n-1);

        c2 += tmp;

        tmp *= n;
        c1 += tmp;

        if( n > 1)
        {
            tmp *= n-1;
            c3 += tmp;

            if( n > 2)
            {
                tmp *= n-1;
                c4 += tmp;
            }
        }


        // field terms
        field_r_mult += c2*cTh*Pn_1-c1*Pn;
        field_z_mult -= c2*Pn_1;

        // how field changes with source direction

        dBdz_grz_rrt_mult += Pn_1*(c1-c2) - cTh*Pn_2*c2;
        dBdz_zzt_mult += c2*(Pn_1+cTh*Pn_2);
        dBdz_zrt_mult += -1.0*c2*Pn_2;
        dBdz_I_mult += -1.0*c2*Pn_1;


        // field gradient
        if( n > 1 )
        {
            grad_I_mult += (c2*cTh*Pn_1-c1*Pn);
            grad_rr_mult += ((c1-c3)*Pn + (2.0*c1-3.0*c2)*cTh*Pn_1 - c2*cThsq*Pn_2);
            grad_zz_mult -= c2*Pn_2;
            grad_rzzr_mult += ((c2-c1)*Pn_1+c2*cTh*Pn_2);
        }

        // how field gradient changes with position
        if( n > 2 )
        {
            gradJacob_zmmz_mult += (c2-c1)*Pn_1+c2*cTh*Pn_2;
            gradJacob_rmmr_mult += (c1-c3)*Pn + (2*c1 - 3*c2)*cTh*Pn_1 - c2*cThsq*Pn_2;
            gradJacob_zrrz_zm_mult += (2.0*c2-c1)*Pn_2+c2*cTh*Pn_3;
            gradJacob_zrrz_rm_mult += (3.0*c1-3.0*c2-c3)*Pn_1 + (2.0*c1-5.0*c2)*cTh*Pn_2 - c2*cThsq*Pn_3;
            gradJacob_I_rm_mult += (c1-c3)*Pn+(2.0*c1-3.0*c2)*cTh*Pn_1-c2*cThsq*Pn_2;
            gradJacob_I_zm_mult += (c2-c1)*Pn_1+c2*cTh*Pn_2;
            gradJacob_zz_zm_mult += -c2*Pn_3;
            gradJacob_zz_rm_mult += (2.0*c2-c1)*Pn_2+c2*cTh*Pn_3;
            gradJacob_rr_zm_mult += (3.0*c1-3.0*c2-c3)*Pn_1+(2.0*c1-5.0*c2)*cTh*Pn_2-c2*cThsq*Pn_3;
            gradJacob_rr_rm_mult += (4.0*c3-3.0*c1-c4)*Pn+(15.0*c2-12.0*c1+3.0*c3)*cTh*Pn_1+(9.0*c2-3.0*c1)*cThsq*Pn_2+c2*cTh*cThsq*Pn_3;
        }

        // how field gradient changes with source direction
        if( n > 1  )
        {
            gradZJacob_I_rm_mult += (c2 - c1)*Pn_1 + (c2*cTh)*Pn_2;
            gradZJacob_I_zm_mult += (-c2)*Pn_2;
            gradZJacob_rr_rm_mult += (3*c1 - 3*c2 - c3)*Pn_1 + (2*c1*cTh - 5*c2*cTh)*Pn_2 + (-c2*cThsq)*Pn_3;
            gradZJacob_rr_zm_mult += (2*c2 - c1)*Pn_2 + (c2*cTh)*Pn_3;
            gradZJacob_zz_rm_mult += (c1 - c2)*Pn_1 + (c1*cTh - 3*c2*cTh)*Pn_2 + (-c2*cThsq)*Pn_3;
            gradZJacob_zz_zm_mult += (2*c2)*Pn_2 + (c2*cTh)*Pn_3;
            gradZJacob_rz_rm_mult += (3*c2*cTh - 3*c1*cTh + c3*cTh)*Pn_1 + (5*c2*cThsq - 2*c1*cThsq)*Pn_2 + (c2*cTh*cThsq)*Pn_3;
            gradZJacob_rz_zm_mult += (c1 - c2)*Pn_1 + (c1*cTh - 3*c2*cTh)*Pn_2 + (-c2*cThsq)*Pn_3;
            gradZJacob_zr_rm_mult += (2*c2 - c1)*Pn_2 + (c2*cTh)*Pn_3;
            gradZJacob_zr_zm_mult += (-c2)*Pn_3;
            gradZJacob_mrrm_mult += (c2 - c1)*Pn_1 + (c2*cTh)*Pn_2;
            gradZJacob_mz_mult += (c1*cTh - c2*cTh)*Pn_1 + (-c2*cThsq)*Pn_2;
            gradZJacob_zm_mult += (-c2)*Pn_2;
        }


    }

    if( pMag > 0 )
    {
        for( coeffIT = src.B_Coeff.begin();
             coeffIT !=src.B_Coeff.end();
             coeffIT ++ )
        {
            double n = coeffIT->order;
            double Pn   = LegandrePolynomial(cTh, n, 0); // Legandre Polynomial at cos(theta)
            double Pn_1 = LegandrePolynomial(cTh, n, 1); // First Derivative of Legandre Polynomial at cos(theta)
            double Pn_2 = LegandrePolynomial(cTh, n, 2); // Second Derivative of Legandre Polynomial at cos(theta)
            double Pn_3 = LegandrePolynomial(cTh, n, 3); // Third Derivative of Legandre Polynomial at cos(theta)

            currentState.value += coeffIT->coeff * 1.0/std::pow(pMag,n+1.0)*Pn;


            double c1 = 0;
            double c2 = 0;
            double c3 = 0;
            double c4 = 0;



            double tmp =  coeffIT->coeff/pow(pMag,n+2);
            c2 += tmp;

            tmp *= (n+1);
            c1 += -1.0*tmp;

            tmp *= n+2;
            c3 += tmp;

            tmp *= n+2;
            c4 += -1.0*tmp;

            // field terms
            field_r_mult += c2*cTh*Pn_1-c1*Pn;
            field_z_mult -= c2*Pn_1;

            // how field changes with source direction

            dBdz_grz_rrt_mult += Pn_1*(c1-c2) - cTh*Pn_2*c2;
            dBdz_zzt_mult += c2*(Pn_1+cTh*Pn_2);
            dBdz_zrt_mult += -1.0*c2*Pn_2;
            dBdz_I_mult += -1.0*c2*Pn_1;


            // field spatial gradient
            grad_I_mult += (c2*cTh*Pn_1-c1*Pn);
            grad_rr_mult += ((c1-c3)*Pn + (2.0*c1-3.0*c2)*cTh*Pn_1 - c2*cThsq*Pn_2);
            grad_zz_mult -= c2*Pn_2;
            grad_rzzr_mult += ((c2-c1)*Pn_1+c2*cTh*Pn_2);

            // how the field spatial gradient changes with position
            gradJacob_zmmz_mult += (c2-c1)*Pn_1+c2*cTh*Pn_2;
            gradJacob_rmmr_mult += (c1-c3)*Pn + (2*c1 - 3*c2)*cTh*Pn_1 - c2*cThsq*Pn_2;
            gradJacob_zrrz_zm_mult += (2.0*c2-c1)*Pn_2+c2*cTh*Pn_3;
            gradJacob_zrrz_rm_mult += (3.0*c1-3.0*c2-c3)*Pn_1 + (2.0*c1-5.0*c2)*cTh*Pn_2 - c2*cThsq*Pn_3;
            gradJacob_I_rm_mult += (c1-c3)*Pn+(2.0*c1-3.0*c2)*cTh*Pn_1-c2*cThsq*Pn_2;
            gradJacob_I_zm_mult += (c2-c1)*Pn_1+c2*cTh*Pn_2;
            gradJacob_zz_zm_mult += -c2*Pn_3;
            gradJacob_zz_rm_mult += (2.0*c2-c1)*Pn_2+c2*cTh*Pn_3;
            gradJacob_rr_zm_mult += (3.0*c1-3.0*c2-c3)*Pn_1+(2.0*c1-5.0*c2)*cTh*Pn_2-c2*cThsq*Pn_3;
            gradJacob_rr_rm_mult += (4.0*c3-3.0*c1-c4)*Pn+(15.0*c2-12.0*c1+3.0*c3)*cTh*Pn_1+(9.0*c2-3.0*c1)*cThsq*Pn_2+c2*cTh*cThsq*Pn_3;


            // how the field spatial gradient changes with source direction
            gradZJacob_I_rm_mult += (c2 - c1)*Pn_1 + (c2*cTh)*Pn_2;
            gradZJacob_I_zm_mult += (-c2)*Pn_2;
            gradZJacob_rr_rm_mult += (3*c1 - 3*c2 - c3)*Pn_1 + (2*c1*cTh - 5*c2*cTh)*Pn_2 + (-c2*cThsq)*Pn_3;
            gradZJacob_rr_zm_mult += (2*c2 - c1)*Pn_2 + (c2*cTh)*Pn_3;
            gradZJacob_zz_rm_mult += (c1 - c2)*Pn_1 + (c1*cTh - 3*c2*cTh)*Pn_2 + (-c2*cThsq)*Pn_3;
            gradZJacob_zz_zm_mult += (2*c2)*Pn_2 + (c2*cTh)*Pn_3;
            gradZJacob_rz_rm_mult += (3*c2*cTh - 3*c1*cTh + c3*cTh)*Pn_1 + (5*c2*cThsq - 2*c1*cThsq)*Pn_2 + (c2*cTh*cThsq)*Pn_3;
            gradZJacob_rz_zm_mult += (c1 - c2)*Pn_1 + (c1*cTh - 3*c2*cTh)*Pn_2 + (-c2*cThsq)*Pn_3;
            gradZJacob_zr_rm_mult += (2*c2 - c1)*Pn_2 + (c2*cTh)*Pn_3;
            gradZJacob_zr_zm_mult += (-c2)*Pn_3;
            gradZJacob_mrrm_mult += (c2 - c1)*Pn_1 + (c2*cTh)*Pn_2;
            gradZJacob_mz_mult += (c1*cTh - c2*cTh)*Pn_1 + (-c2*cThsq)*Pn_2;
            gradZJacob_zm_mult += (-c2)*Pn_2;

            //            // The second zHat deritive
            //            zSqJacob_I_rm_mult += (-c2)*Pn_2;
            //            zSqJacob_I_zm_mult += c2*Pn_1 + (c2*cTh)*Pn_2;
            //            zSqJacob_rr_rm_mult += (2*c2 - c1)*Pn_2 + (c2*cTh)*Pn_3;
            //            zSqJacob_rr_zm_mult += (c1 - c2)*Pn_1 + (c1*cTh - 3*c2*cTh)*Pn_2 + (-c2*cThsq)*Pn_3;
            //            zSqJacob_zz_rm_mult += c2*Pn_2 + (c2*cTh)*Pn_3;
            //            zSqJacob_zz_zm_mult += (-2*c2)*Pn_1 + (-4*c2*cTh)*Pn_2 + (-c2*cThsq)*Pn_3;
            //            zSqJacob_rz_rm_mult += (c1*cTh - 2*c2*cTh)*Pn_2 + (-c2*cThsq)*Pn_3;
            //            zSqJacob_rz_zm_mult += (2*c2*cTh - 2*c1*cTh)*Pn_1 + (4*c2*cThsq - c1*cThsq)*Pn_2 + (c2*cThsq*cTh)*Pn_3;
            //            zSqJacob_zr_rm_mult += (-c2)*Pn_3;
            //            zSqJacob_zr_zm_mult += (2*c2)*Pn_2 + (c2*cTh)*Pn_3;
            //            zSqJacob_mr_mult += (-c2)*Pn_2;
            //            zSqJacob_rm_mult += (c1*cTh - c2*cTh)*Pn_1 + (-c2*cThsq)*Pn_2;
            //            zSqJacob_mz_mult += (c2*cTh)*Pn_2;
            //            zSqJacob_zm_mult += c2*Pn_1 + (c2*cTh)*Pn_2;

        }


    }



    // calculate field
    currentState.firstSpatialDerivative    += field_r_mult*p + field_z_mult*z;

    // claculate field spatial gradient
    Eigen::Matrix3d GradMat =  (grad_I_mult/pMag)*I
            +(grad_rr_mult/pMag)*(ppt)
            +(grad_zz_mult/pMag)*(zzt)
            +(grad_rzzr_mult/pMag)*(pzt_zpt);

    currentState.secondSpatialDerivative +=GradMat;
    currentState.firstSpatialDerivative_SourcePositionDerivative.push_back(-GradMat);



    // calculate how field changes with source direction

    currentState.firstSpatialDerivative_SourceHeadingDerivative.push_back( dBdz_grz_rrt_mult*(cTh*rzt-ppt)
                                                                           + dBdz_zzt_mult*zzt
                                                                           + dBdz_zrt_mult*zrt
                                                                           + dBdz_I_mult *I );



    // initialze terms for 3rd Derivatives
    Eigen::Matrix3d rmt = p*X.transpose();
    Eigen::Matrix3d mrt = X*p.transpose();
    Eigen::Matrix3d rm_mr = rmt + mrt;

    Eigen::Matrix3d mzt = X*z.transpose();
    Eigen::Matrix3d zmt = z*X.transpose();
    Eigen::Matrix3d zm_mz = zmt+mzt;

    double r_dot_m = p(0);
    double z_dot_m = z(0);

    Eigen::Matrix<double,5,3> gradientSpatialDer;
    Eigen::Matrix<double,5,3> gradientHeadingDer;
    Eigen::Matrix<double,9,3> fieldHedingSecDer; fieldHedingSecDer.setZero(9,3);

    // calculate how the field spatial Derivative changes with position (first 3 rows)
    gradientSpatialDer.topRows<3>() = zm_mz*(gradJacob_zmmz_mult/pMagSq)
            +rm_mr*(gradJacob_rmmr_mult/pMagSq)
            +pzt_zpt*((gradJacob_zrrz_zm_mult*z_dot_m+gradJacob_zrrz_rm_mult*r_dot_m)/pMagSq)
            +I*((gradJacob_I_zm_mult*z_dot_m+gradJacob_I_rm_mult*r_dot_m)/pMagSq)
            +zzt*((gradJacob_zz_rm_mult*r_dot_m + gradJacob_zz_zm_mult*z_dot_m)/pMagSq)
            +ppt*((gradJacob_rr_rm_mult*r_dot_m + gradJacob_rr_zm_mult*z_dot_m)/pMagSq);


    // calculate how the field spatial Derivative changes with source direction (first 3 rows)
    gradientHeadingDer.topRows<3>() = I*((gradZJacob_I_rm_mult*r_dot_m+gradZJacob_I_zm_mult*z_dot_m)/pMag)
            +ppt*((gradZJacob_rr_rm_mult*r_dot_m+gradZJacob_rr_zm_mult*z_dot_m)/pMag)
            +zzt*((gradZJacob_zz_rm_mult*r_dot_m+gradZJacob_zz_zm_mult*z_dot_m)/pMag)
            +rzt*((gradZJacob_rz_rm_mult*r_dot_m+gradZJacob_rz_zm_mult*z_dot_m)/pMag)
            +zrt*((gradZJacob_zr_rm_mult*r_dot_m+gradZJacob_zr_zm_mult*z_dot_m)/pMag)
            +(rm_mr)*((gradZJacob_mrrm_mult)/pMag)
            +zmt*((gradZJacob_zm_mult)/pMag)
            +mzt*((gradZJacob_mz_mult)/pMag);

    //    fieldHedingSecDer.topRows<3>() = I*((zSqJacob_I_rm_mult*r_dot_m+zSqJacob_I_zm_mult*z_dot_m))
    //            +ppt*((zSqJacob_rr_rm_mult*r_dot_m+zSqJacob_rr_zm_mult*z_dot_m))
    //            +zzt*((zSqJacob_zz_rm_mult*r_dot_m+zSqJacob_zz_zm_mult*z_dot_m))
    //            +rzt*((zSqJacob_rz_rm_mult*r_dot_m+zSqJacob_rz_zm_mult*z_dot_m))
    //            +zrt*((zSqJacob_zr_rm_mult*r_dot_m+zSqJacob_zr_zm_mult*z_dot_m))
    //            +mrt*(zSqJacob_mr_mult) + rmt*(zSqJacob_rm_mult)
    //            +zmt*((zSqJacob_zm_mult))
    //            +mzt*((zSqJacob_mz_mult));



    // Y Component Time
    rmt = p*Y.transpose();
    mrt = Y*p.transpose();
    rm_mr = rmt+mrt;
    mzt = Y*z.transpose();
    zmt = z*Y.transpose();
    zm_mz = zmt+mzt;

    r_dot_m = p(1);
    z_dot_m = z(1);

    // calculate how the field spatial Derivative changes with position (last 2 rows)
    gradientSpatialDer.bottomRows<2>() = zm_mz.bottomRows<2>()*(gradJacob_zmmz_mult/pMagSq)
            +rm_mr.bottomRows<2>()*(gradJacob_rmmr_mult/pMagSq)
            +pzt_zpt.bottomRows<2>()*((gradJacob_zrrz_zm_mult*z_dot_m+gradJacob_zrrz_rm_mult*r_dot_m)/pMagSq)
            +I.bottomRows<2>()*((gradJacob_I_zm_mult*z_dot_m+gradJacob_I_rm_mult*r_dot_m)/pMagSq)
            +zzt.bottomRows<2>()*((gradJacob_zz_rm_mult*r_dot_m + gradJacob_zz_zm_mult*z_dot_m)/pMagSq)
            +ppt.bottomRows<2>()*((gradJacob_rr_rm_mult*r_dot_m + gradJacob_rr_zm_mult*z_dot_m)/pMagSq);


    // calculate how the field spatial Derivative changes with source direction (last 2 rows)
    gradientHeadingDer.bottomRows<2>() = I.bottomRows<2>()*((gradZJacob_I_rm_mult* r_dot_m+gradZJacob_I_zm_mult* z_dot_m)/pMag)
            +ppt.bottomRows<2>()*((gradZJacob_rr_rm_mult*r_dot_m+gradZJacob_rr_zm_mult*z_dot_m)/pMag)
            +zzt.bottomRows<2>()*((gradZJacob_zz_rm_mult*r_dot_m+gradZJacob_zz_zm_mult*z_dot_m)/pMag)
            +rzt.bottomRows<2>()*((gradZJacob_rz_rm_mult*r_dot_m+gradZJacob_rz_zm_mult*z_dot_m)/pMag)
            +zrt.bottomRows<2>()*((gradZJacob_zr_rm_mult*r_dot_m+gradZJacob_zr_zm_mult*z_dot_m)/pMag)
            +(rm_mr).bottomRows<2>()*((gradZJacob_mrrm_mult)/pMag)
            +zmt.bottomRows<2>()*((gradZJacob_zm_mult)/pMag)
            +mzt.bottomRows<2>()*((gradZJacob_mz_mult)/pMag);

    //    fieldHedingSecDer.block<3,3>(3,0) = I*((zSqJacob_I_rm_mult*r_dot_m+zSqJacob_I_zm_mult*z_dot_m))
    //            +ppt*((zSqJacob_rr_rm_mult*r_dot_m+zSqJacob_rr_zm_mult*z_dot_m))
    //            +zzt*((zSqJacob_zz_rm_mult*r_dot_m+zSqJacob_zz_zm_mult*z_dot_m))
    //            +rzt*((zSqJacob_rz_rm_mult*r_dot_m+zSqJacob_rz_zm_mult*z_dot_m))
    //            +zrt*((zSqJacob_zr_rm_mult*r_dot_m+zSqJacob_zr_zm_mult*z_dot_m))
    //            +mrt*(zSqJacob_mr_mult) + rmt*(zSqJacob_rm_mult)
    //            +zmt*((zSqJacob_zm_mult))
    //            +mzt*((zSqJacob_mz_mult));


    //    Eigen::Vector3d Z(0,0,1);
    //    rmt = p*Z.transpose();
    //    mrt = Z*p.transpose();
    //    zmt = z*Z.transpose();
    //    mzt = Z*z.transpose();
    //    r_dot_m = p(2);
    //    z_dot_m = z(2);

    //    fieldHedingSecDer.bottomRows<3>() = I*((zSqJacob_I_rm_mult*r_dot_m+zSqJacob_I_zm_mult*z_dot_m))
    //            +ppt*((zSqJacob_rr_rm_mult*r_dot_m+zSqJacob_rr_zm_mult*z_dot_m))
    //            +zzt*((zSqJacob_zz_rm_mult*r_dot_m+zSqJacob_zz_zm_mult*z_dot_m))
    //            +rzt*((zSqJacob_rz_rm_mult*r_dot_m+zSqJacob_rz_zm_mult*z_dot_m))
    //            +zrt*((zSqJacob_zr_rm_mult*r_dot_m+zSqJacob_zr_zm_mult*z_dot_m))
    //            +mrt*(zSqJacob_mr_mult) + rmt*(zSqJacob_rm_mult)
    //            +zmt*((zSqJacob_zm_mult))
    //            +mzt*((zSqJacob_mz_mult));



    currentState.thirdSpatialDerivative += gradientSpatialDer;


    currentState.secondSpatialDerivative_SourcePositionDerivative.push_back(gradientSpatialDer);
    currentState.secondSpatialDerivative_SourceHeadingDerivative.push_back(-gradientHeadingDer);
    //    currentState.firstSpatialDerivative_secondSourceHeadingDerivative.push_back(fieldHedingSecDer);



    return ;
}

// This fuction returns the field coeffient Derivatives and direction Derivatives assuming a current of 1 at the position specified in the work space
ScalorPotentialCalibrationJacobians ScalorPotential::srcCalibrationInformation(const Eigen::Vector3d& position, unsigned int srcNum ) const
{
    assert( srcNum < srcList.size() );

    ScalorPotentialCalibrationJacobians currentStateJacob;

    // Use forward function to get data for B, dBdps, and dBdzs
    ScalorPotentialState currentState;
    srcFieldGradient(position, srcList[srcNum], currentState );
    currentStateJacob.firstSpatialDerivative = currentState.firstSpatialDerivative;
    currentStateJacob.firstSpatialDerivative_SourcePositionDerivative = currentState.firstSpatialDerivative_SourcePositionDerivative.front();
    currentStateJacob.firstSpatialDerivative_SourceHeadingDerivative = currentState.firstSpatialDerivative_SourceHeadingDerivative.front();
    currentStateJacob.secondSpatialDerivative_SourceHeadingDerivative = currentState.secondSpatialDerivative_SourceHeadingDerivative.front();
    currentStateJacob.secondSpatialDerivative_SourcePositionDerivative = currentState.secondSpatialDerivative_SourcePositionDerivative.front();

    // Zero out other coefficients
    currentStateJacob.firstSpatialDerivative_A_CoeffDerivative.setZero(3,srcList[srcNum].A_Coeff.size());
    currentStateJacob.firstSpatialDerivative_B_CoeffDerivative.setZero(3,srcList[srcNum].B_Coeff.size());
    currentStateJacob.firstSpatialDerivative_dA_dHeading.setZero(9,srcList[srcNum].A_Coeff.size());
    currentStateJacob.firstSpatialDerivative_dB_dHeading.setZero(9,srcList[srcNum].B_Coeff.size());
    currentStateJacob.firstSpatialDerivative_dA_dPosition.setZero(9,srcList[srcNum].A_Coeff.size());
    currentStateJacob.firstSpatialDerivative_dB_dPosition.setZero(9,srcList[srcNum].B_Coeff.size());


    Eigen::Vector3d p(position - srcList[srcNum].srcPosition);
    Eigen::Vector3d z(srcList[srcNum].srcDirection.normalized());
    double pMag = p.norm();
    if( pMag > 0 )
    {
        p/=pMag;
    }

    // define Cosine Theta
    double cTh = p.dot(z);
    double cThsq = cTh*cTh;

    Eigen::Matrix3d pzt_zpt( p*z.transpose()+z*p.transpose() );
    Eigen::Matrix3d ppt( p*p.transpose() );// = p*p.transpose(); // p*p'
    Eigen::Matrix3d zzt(z *z.transpose() );// = src.srcDirection *src.srcDirection.transpose(); // z*z'
    Eigen::Matrix3d rzt = p*z.transpose();
    Eigen::Matrix3d zrt = z*p.transpose();

    // some constant vectors and matricies
    Eigen::Matrix3d I(Eigen::Matrix3d::Identity());


    std::vector<srcCoeff>::const_iterator coeffIT;
    int colNum = 0;

    // The second zHat deritive
    double zSqJacob_I_rm_mult = 0;
    double zSqJacob_I_zm_mult = 0;
    double zSqJacob_rr_rm_mult = 0;
    double zSqJacob_rr_zm_mult = 0;
    double zSqJacob_zz_rm_mult = 0;
    double zSqJacob_zz_zm_mult = 0;
    double zSqJacob_rz_rm_mult = 0;
    double zSqJacob_rz_zm_mult = 0;
    double zSqJacob_zr_rm_mult = 0;
    double zSqJacob_zr_zm_mult = 0;
    double zSqJacob_mr_mult = 0;
    double zSqJacob_rm_mult = 0;
    double zSqJacob_mz_mult = 0;
    double zSqJacob_zm_mult = 0;

    for( coeffIT = srcList[srcNum].A_Coeff.begin();
         coeffIT !=srcList[srcNum].A_Coeff.end();
         coeffIT ++ )
    {
        double n = coeffIT->order;

        currentState.value += coeffIT->coeff * std::pow(pMag,n);

        if( n==0 )
        {
            colNum ++;
            continue;
        }

        double Pn   = LegandrePolynomial(cTh, n, 0); // Legandre Polynomial at cos(theta)
        double Pn_1 = LegandrePolynomial(cTh, n, 1); // First Derivative of Legandre Polynomial at cos(theta)
        double Pn_2 = LegandrePolynomial(cTh, n, 2); // Second Derivative of Legandre Polynomial at cos(theta)
        double Pn_3 = LegandrePolynomial(cTh, n, 3); // Second Derivative of Legandre Polynomial at cos(theta)


        double c1 = 0;
        double c2 = 0;
        double c3 = 0;

        double tmp =  1;


        if( n > 1 )
            tmp *= pow(pMag,n-1);

        c2 += tmp;

        tmp *= n;
        c1 += tmp;

        if( n > 1)
        {
            tmp *= n-1;
            c3 += tmp;
        }


        // field terms

        double field_r_mult = c2*cTh*Pn_1-c1*Pn;
        double field_z_mult = -c2*Pn_1;

        currentStateJacob.firstSpatialDerivative_A_CoeffDerivative.col(colNum) = field_r_mult*p + field_z_mult*z;


        // how field changes with source direction
        double dBdz_grz_rrt_mult = Pn_1*(c1-c2) - cTh*Pn_2*c2;
        double dBdz_zzt_mult = c2*(Pn_1+cTh*Pn_2);
        double dBdz_zrt_mult = -1.0*c2*Pn_2;
        double dBdz_I_mult = -1.0*c2*Pn_1;




        // field gradient
        double grad_I_mult = (c2*cTh*Pn_1-c1*Pn);
        double grad_rr_mult = ((c1-c3)*Pn + (2.0*c1-3.0*c2)*cTh*Pn_1 - c2*cThsq*Pn_2);
        double grad_zz_mult = -c2*Pn_2;
        double grad_rzzr_mult = ((c2-c1)*Pn_1+c2*cTh*Pn_2);

        // claculate field spatial gradient
        Eigen::Matrix3d GradMat =  -1*((grad_I_mult/pMag)*I
                                       +(grad_rr_mult/pMag)*(ppt)
                                       +(grad_zz_mult/pMag)*(zzt)
                                       +(grad_rzzr_mult/pMag)*(pzt_zpt));



        Eigen::Matrix3d HeadMat = dBdz_grz_rrt_mult*(cTh*rzt-ppt)
                + dBdz_zzt_mult*zzt
                + dBdz_zrt_mult*zrt
                + dBdz_I_mult *I;

        // this is needed because we actually want d(B'X)/dAdHeading

        //HeadMat.transposeInPlace();


        for( int i=0; i<3; i++ )
        {
            currentStateJacob.firstSpatialDerivative_dA_dPosition.block<3,1>(i*3,colNum) = GradMat.transpose().col(i);
            currentStateJacob.firstSpatialDerivative_dA_dHeading.block<3,1>(i*3,colNum) = HeadMat.col(i);
        }

        // The second zHat deritive
        c1 *= coeffIT->coeff;
        c2 *= coeffIT->coeff;
        c3 *= coeffIT->coeff;
        zSqJacob_I_rm_mult += ((-c2)*Pn_2);
        zSqJacob_I_zm_mult += c2*Pn_1 + (c2*cTh)*Pn_2;
        zSqJacob_rr_rm_mult += (2*c2 - c1)*Pn_2 + (c2*cTh)*Pn_3;
        zSqJacob_rr_zm_mult += (c1 - c2)*Pn_1 + (c1*cTh - 3*c2*cTh)*Pn_2 + (-c2*cThsq)*Pn_3;
        zSqJacob_zz_rm_mult += c2*Pn_2 + (c2*cTh)*Pn_3;
        zSqJacob_zz_zm_mult += (-2*c2)*Pn_1 + (-4*c2*cTh)*Pn_2 + (-c2*cThsq)*Pn_3;
        zSqJacob_rz_rm_mult += (c1*cTh - 2*c2*cTh)*Pn_2 + (-c2*cThsq)*Pn_3;
        zSqJacob_rz_zm_mult += (2*c2*cTh - 2*c1*cTh)*Pn_1 + (4*c2*cThsq - c1*cThsq)*Pn_2 + (c2*cThsq*cTh)*Pn_3;
        zSqJacob_zr_rm_mult += (-c2)*Pn_3;
        zSqJacob_zr_zm_mult += (2*c2)*Pn_2 + (c2*cTh)*Pn_3;
        zSqJacob_mr_mult += (-c2)*Pn_2;
        zSqJacob_rm_mult += (c1*cTh - c2*cTh)*Pn_1 + (-c2*cThsq)*Pn_2;
        zSqJacob_mz_mult += (c2*cTh)*Pn_2;
        zSqJacob_zm_mult += c2*Pn_1 + (c2*cTh)*Pn_2;

        colNum++;


    }

    if( pMag > 0 )
    {
        colNum = 0;
        for( coeffIT = srcList[srcNum].B_Coeff.begin();
             coeffIT !=srcList[srcNum].B_Coeff.end();
             coeffIT ++ )
        {
            double n = coeffIT->order;

            currentState.value += coeffIT->coeff * 1.0/std::pow(pMag,n+1.0);

            double Pn   = LegandrePolynomial(cTh, n, 0); // Legandre Polynomial at cos(theta)
            double Pn_1 = LegandrePolynomial(cTh, n, 1); // First Derivative of Legandre Polynomial at cos(theta)
            double Pn_2 = LegandrePolynomial(cTh, n, 2); // Second Derivative of Legandre Polynomial at cos(theta)
            double Pn_3 = LegandrePolynomial(cTh, n, 3); // Second Derivative of Legandre Polynomial at cos(theta)


            double c1 = 0;
            double c2 = 0;
            double c3 = 0;

            double tmp =  1.0/pow(pMag,n+2);
            c2 += tmp;

            tmp *= (n+1);
            c1 += -1.0*tmp;

            tmp *= n+2;
            c3 += tmp;


            // field terms
            double field_r_mult = c2*cTh*Pn_1-c1*Pn;
            double field_z_mult = -c2*Pn_1;

            currentStateJacob.firstSpatialDerivative_B_CoeffDerivative.col(colNum) = field_r_mult*p + field_z_mult*z;


            // how field changes with source direction
            double dBdz_grz_rrt_mult = Pn_1*(c1-c2) - cTh*Pn_2*c2;
            double dBdz_zzt_mult = c2*(Pn_1+cTh*Pn_2);
            double dBdz_zrt_mult = -1.0*c2*Pn_2;
            double dBdz_I_mult = -1.0*c2*Pn_1;




            // field gradient
            double grad_I_mult = (c2*cTh*Pn_1-c1*Pn);
            double grad_rr_mult = ((c1-c3)*Pn + (2.0*c1-3.0*c2)*cTh*Pn_1 - c2*cThsq*Pn_2);
            double grad_zz_mult = -c2*Pn_2;
            double grad_rzzr_mult = ((c2-c1)*Pn_1+c2*cTh*Pn_2);

            // claculate field spatial gradient
            Eigen::Matrix3d GradMat =  -1*((grad_I_mult/pMag)*I
                                           +(grad_rr_mult/pMag)*(ppt)
                                           +(grad_zz_mult/pMag)*(zzt)
                                           +(grad_rzzr_mult/pMag)*(pzt_zpt));



            Eigen::Matrix3d HeadMat = dBdz_grz_rrt_mult*(cTh*rzt-ppt)
                    + dBdz_zzt_mult*zzt
                    + dBdz_zrt_mult*zrt
                    + dBdz_I_mult *I;

            // this is needed because we actually want d(B'X)/dBdHeading
            //HeadMat.transposeInPlace();

            for( int i=0; i<3; i++ )
            {
                currentStateJacob.firstSpatialDerivative_dB_dPosition.block<3,1>(i*3,colNum) = GradMat.col(i);
                currentStateJacob.firstSpatialDerivative_dB_dHeading.block<3,1>(i*3,colNum) = HeadMat.col(i);
            }

            c1 *= coeffIT->coeff;
            c2 *= coeffIT->coeff;
            c3 *= coeffIT->coeff;
            zSqJacob_I_rm_mult += ((-c2)*Pn_2);
            zSqJacob_I_zm_mult += c2*Pn_1 + (c2*cTh)*Pn_2;
            zSqJacob_rr_rm_mult += (2*c2 - c1)*Pn_2 + (c2*cTh)*Pn_3;
            zSqJacob_rr_zm_mult += (c1 - c2)*Pn_1 + (c1*cTh - 3*c2*cTh)*Pn_2 + (-c2*cThsq)*Pn_3;
            zSqJacob_zz_rm_mult += c2*Pn_2 + (c2*cTh)*Pn_3;
            zSqJacob_zz_zm_mult += (-2*c2)*Pn_1 + (-4*c2*cTh)*Pn_2 + (-c2*cThsq)*Pn_3;
            zSqJacob_rz_rm_mult += (c1*cTh - 2*c2*cTh)*Pn_2 + (-c2*cThsq)*Pn_3;
            zSqJacob_rz_zm_mult += (2*c2*cTh - 2*c1*cTh)*Pn_1 + (4*c2*cThsq - c1*cThsq)*Pn_2 + (c2*cThsq*cTh)*Pn_3;
            zSqJacob_zr_rm_mult += (-c2)*Pn_3;
            zSqJacob_zr_zm_mult += (2*c2)*Pn_2 + (c2*cTh)*Pn_3;
            zSqJacob_mr_mult += (-c2)*Pn_2;
            zSqJacob_rm_mult += (c1*cTh - c2*cTh)*Pn_1 + (-c2*cThsq)*Pn_2;
            zSqJacob_mz_mult += (c2*cTh)*Pn_2;
            zSqJacob_zm_mult += c2*Pn_1 + (c2*cTh)*Pn_2;



            colNum++;

        }


    }

    Eigen::Matrix<double,9,3> fieldHedingSecDer;
    double r_dot_m = p(0);
    double z_dot_m = z(0);
    Eigen::Vector3d X(1,0,0),Y(0,1,0),Z(0,0,1);
    Eigen::Matrix3d  rmt = p*X.transpose();
    Eigen::Matrix3d  mrt = X*p.transpose();
    Eigen::Matrix3d  mzt = X*z.transpose();
    Eigen::Matrix3d  zmt = z*X.transpose();

    fieldHedingSecDer.topRows<3>() = I*((zSqJacob_I_rm_mult*r_dot_m+zSqJacob_I_zm_mult*z_dot_m))
            +ppt*((zSqJacob_rr_rm_mult*r_dot_m+zSqJacob_rr_zm_mult*z_dot_m))
            +zzt*((zSqJacob_zz_rm_mult*r_dot_m+zSqJacob_zz_zm_mult*z_dot_m))
            +rzt*((zSqJacob_rz_rm_mult*r_dot_m+zSqJacob_rz_zm_mult*z_dot_m))
            +zrt*((zSqJacob_zr_rm_mult*r_dot_m+zSqJacob_zr_zm_mult*z_dot_m))
            +mrt*(zSqJacob_mr_mult) + rmt*(zSqJacob_rm_mult)
            +zmt*((zSqJacob_zm_mult))
            +mzt*((zSqJacob_mz_mult));



    // Y Component Time
    rmt = p*Y.transpose();
    mrt = Y*p.transpose();
    mzt = Y*z.transpose();
    zmt = z*Y.transpose();
    r_dot_m = p(1);
    z_dot_m = z(1);


    fieldHedingSecDer.block<3,3>(3,0) = I*((zSqJacob_I_rm_mult*r_dot_m+zSqJacob_I_zm_mult*z_dot_m))
            +ppt*((zSqJacob_rr_rm_mult*r_dot_m+zSqJacob_rr_zm_mult*z_dot_m))
            +zzt*((zSqJacob_zz_rm_mult*r_dot_m+zSqJacob_zz_zm_mult*z_dot_m))
            +rzt*((zSqJacob_rz_rm_mult*r_dot_m+zSqJacob_rz_zm_mult*z_dot_m))
            +zrt*((zSqJacob_zr_rm_mult*r_dot_m+zSqJacob_zr_zm_mult*z_dot_m))
            +mrt*(zSqJacob_mr_mult) + rmt*(zSqJacob_rm_mult)
            +zmt*((zSqJacob_zm_mult))
            +mzt*((zSqJacob_mz_mult));


    rmt = p*Z.transpose();
    mrt = Z*p.transpose();
    zmt = z*Z.transpose();
    mzt = Z*z.transpose();
    r_dot_m = p(2);
    z_dot_m = z(2);

    fieldHedingSecDer.bottomRows<3>() = I*((zSqJacob_I_rm_mult*r_dot_m+zSqJacob_I_zm_mult*z_dot_m))
            +ppt*((zSqJacob_rr_rm_mult*r_dot_m+zSqJacob_rr_zm_mult*z_dot_m))
            +zzt*((zSqJacob_zz_rm_mult*r_dot_m+zSqJacob_zz_zm_mult*z_dot_m))
            +rzt*((zSqJacob_rz_rm_mult*r_dot_m+zSqJacob_rz_zm_mult*z_dot_m))
            +zrt*((zSqJacob_zr_rm_mult*r_dot_m+zSqJacob_zr_zm_mult*z_dot_m))
            +mrt*(zSqJacob_mr_mult) + rmt*(zSqJacob_rm_mult)
            +zmt*((zSqJacob_zm_mult))
            +mzt*((zSqJacob_mz_mult));

    currentStateJacob.firstSpatialDerivative_secondSourceHeadingDerivative = fieldHedingSecDer;

    return currentStateJacob;
}



// This function returns the value of the legandre polynomial of order "order" with Derivative der at position x in [-1,1]
double ScalorPotential::LegandrePolynomial(double x, int order, int der)
{


    if( !(x >= -1.0 && x <=1.0) )
    {
        //cout << "ERROR: LegandrePolynomial Call with X outside of [-1,1] inclusive range.  Value is: " << setprecision(12)<< x << " Clamping to range..." << endl;
        x = (x>=0)?std::min(x,1.0):std::max(x,-1.0);

    }
    if( std::isnan(x) || std::isinf(x) )
        x = 0;



    //assert( x >= -1 && x <=1 ); // polynomial must be between -1 and 1 to be valid

    double value = 0;
    if( order >= der )
    {
        for( double k = der; k <= order; k++ )
        {
            double multiplier = 1;
            // N choose K
            for( double j=1.0; j<= k; j++ )
            {
                multiplier *= (order+1.0-j)/j;
            }
            // (n+k-1)/2.0 choose n  / 2^n
            for( double j=1.0; j<=order; j++)
            {
                multiplier *= (order+k-2.0*j+1)/j;
            }
            // Derivative coeff product
            for(double j=0.0;j<der;j++)
            {
                multiplier *= (k-j);
            }
            value += pow(x,k-der)*multiplier;
        }
    }
    return value;
}


// **********************   CALIBRATION FUNCTIONS BELOW ************************** //
//void ScalorPotential::packCalibrationState(Eigen::VectorXd& stateVector) const
//{

//    stateVector.setZero(getNumCalibrationParameters(),1);
//    int vectorIndex = 0;

//    std::vector<srcStruct>::const_iterator srcIT = srcList.begin();

//    for( ; srcIT != srcList.end(); srcIT ++ )
//    {
//        std::vector<srcCoeff>::const_iterator coeffIT = srcIT->A_Coeff.begin();

//        // unpack A coeff
//        for( ; coeffIT != srcIT->A_Coeff.end(); coeffIT ++, vectorIndex ++ )
//        {
//            stateVector(vectorIndex) = coeffIT->coeff;
//        }

//        // unpack B coeff
//        for( srcIT->B_Coeff.begin(); coeffIT != srcIT->A_Coeff.end(); coeffIT ++, vectorIndex ++ )
//        {
//            stateVector(vectorIndex) = coeffIT->coeff;
//        }

//        // unpack direction
//        for( int i=0; i<3; i++, vectorIndex ++ )
//        {
//            stateVector(vectorIndex) = srcIT->srcDirection(i);
//        }

//        // unpack position
//        for( int i=0; i<3; i++, vectorIndex ++ )
//        {
//            stateVector(vectorIndex) = srcIT->srcPosition(i);
//        }
//    }
//}
//void ScalorPotential::unpackCalibrationState(Eigen::VectorXd& stateVector)
//{
//    assert( ('State vector provided for unpacking is the worng size', stateVector.rows() == getNumCalibrationParameters() ));

//    int vectorIndex = 0;

//    std::vector<srcStruct>::iterator srcIT = srcList.begin();

//    for( ; srcIT != srcList.end(); srcIT ++ )
//    {
//        std::vector<srcCoeff>::iterator coeffIT = srcIT->A_Coeff.begin();

//        // unpack A coeff
//        for( ; coeffIT != srcIT->A_Coeff.end(); coeffIT ++, vectorIndex ++ )
//        {
//            coeffIT->coeff = stateVector(vectorIndex);
//        }

//        // unpack B coeff
//        for( srcIT->B_Coeff.begin(); coeffIT != srcIT->A_Coeff.end(); coeffIT ++, vectorIndex ++ )
//        {
//            coeffIT->coeff = stateVector(vectorIndex);
//        }

//        // unpack direction
//        for( int i=0; i<3; i++, vectorIndex ++ )
//        {
//            srcIT->srcDirection(i) = stateVector(vectorIndex);
//        }
//        assert(('Direction has zero length!',srcIT->srcDirection.norm() != 0 ));
//        srcIT->srcDirection.normalize();

//        // unpack position
//        for( int i=0; i<3; i++, vectorIndex ++ )
//        {
//            srcIT->srcPosition(i) = stateVector(vectorIndex);
//        }
//    }

//    assert( vectorIndex == getNumCalibrationParameters());
//}

int ScalorPotential::getNumCalibrationParameters( int srcNum ) const
{
    if( srcNum == -1 )
        return numCalParameters;
    else
        return srcList[srcNum].A_Coeff.size() + srcList[srcNum].B_Coeff.size() +6;
}

