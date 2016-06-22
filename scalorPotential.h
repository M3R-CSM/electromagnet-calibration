#ifndef ScalorPotential_H
#define ScalorPotential_H
/*
 ScalorPotential is a calss that returns the spherical scalor potential information at a point.

    Public Functions:
        getState( position, sourceNumber ) :  Returns the ScalorPotentialState at the specified point. \see ScalorPotentialState

        getValue( position, sourceNumber ) : Returns the value of the scalor potential

        getGradient( position, sourceNumber ) : Returns the first spatial deritive of the scalor potential

        getNumberOfSources( ): returns the current number of sources

        getSourceStruct( sourceNumber ): returns the desired source struct

        setSourceStruct( sourceNumber, newSrc): sets the source desired to a new value

        removeSourceStruct( sourceNumber ); removes the source from the list

        remapSecondDerivativeVec( derivative vector ) : maps the 5x1 vector to a 3x3 symmetric zero-trace spatial gradient matrix

        remapSecondDerivativeMat( dervative Matrix) : maps the 3x3 symmetric zero-trace spacial gradient matrix to the 5x1 vector

        ScalorPotential(): Initialzie as empty

        ScalorPotential( const std::vector<srcStruct>& srcList ): Initialize with a source list

    Protected Functions:
        LegandrePolynomial: Calculates the Legandre Polynomial of a given order deritivative X <= (-1,1)
        srcFieldGradient: updates the current potential state with a particular source's contribution
        srcCalibrationInformation: generates the matricies required for calibration for a given source

        packCalibrationState: packs the calibration state vector
        unpackCalibrationState: unpacks the calibration state vector
        getNumCalibrationParameters:  returns the calibration matrix size for this source
        packCalibrationErrorVector:  packs the calibration error vector for this source
        packCalibrationMatrix:  packs the calibration matrix for this source


    Libraries:
    This class requires the Eigen library version 3 or higher. http://eigen.tuxfamily.org/
        Make sure to update the #include <Eigen/Core> line for your library specifics

    Files:
        This class requires the following files to compile
            ScalorPotential.h
            ScalorPotential.cpp

**/

#include <Eigen/Core>

#include <limits>
#include <vector>
#include <assert.h>

class ElectromagnetCalibration;

typedef Eigen::Matrix<double,5,1> Vector5d;

struct ScalorPotentialState
{
    double value; ///< The value of the scalor potential
    Eigen::Vector3d firstSpatialDerivative; ///< Field, the gradient of the potential.
    Eigen::Matrix<double,3,3> secondSpatialDerivative; ///< Field spatial gradient
    Eigen::Matrix<double,5,3> thirdSpatialDerivative;  ///< How the field spatial gradient changes with position, assumes potentials of order > 1. (no sources or sinks)

    std::vector<Eigen::MatrixXd> firstSpatialDerivative_SourceHeadingDerivative;  ///< A list of 3x3 matricies describing how the field changes with the source heading
    std::vector<Eigen::MatrixXd> firstSpatialDerivative_SourcePositionDerivative; ///< A list of 3x3 matricies describing how the field changes with the source position

    std::vector<Eigen::MatrixXd> secondSpatialDerivative_SourcePositionDerivative; ///< A list of 5x3 matricies describing how the field spatial gradient changes with the source position
    std::vector<Eigen::MatrixXd> secondSpatialDerivative_SourceHeadingDerivative;  ///< A list of 5x3 matricies describing how the field spatial gradient changes with the source heading


    ScalorPotentialState();

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

struct ScalorPotentialCalibrationJacobians
{
    Eigen::Vector3d firstSpatialDerivative; ///< Field, the gradient of the potential.

    Eigen::MatrixXd firstSpatialDerivative_A_CoeffDerivative; ///< How the field changes with the A coefficients
    Eigen::MatrixXd firstSpatialDerivative_B_CoeffDerivative; ///< How the field changes with the B coefficients
    Eigen::MatrixXd firstSpatialDerivative_dA_dHeading; ///< The second deritive of the field with respect to the A coefficients and with Heading
    Eigen::MatrixXd firstSpatialDerivative_dB_dHeading; ///< The second deritive of the field with respect to the B coefficients and with Heading
    Eigen::MatrixXd firstSpatialDerivative_dA_dPosition; ///< The second deritive of the field with respect to the A coefficients and with Position
    Eigen::MatrixXd firstSpatialDerivative_dB_dPosition; ///< The second deritive of the field with respect to the B coefficients and with Position

    Eigen::Matrix3d firstSpatialDerivative_SourcePositionDerivative; ///< Field spatial gradient
    Eigen::Matrix3d firstSpatialDerivative_SourceHeadingDerivative; ///< How the field changes with the source heading

    Eigen::Matrix<double, 5,3 > secondSpatialDerivative_SourcePositionDerivative; ///< 5x3 matrix describing how the field spatial gradient changes with the source position
    Eigen::Matrix<double, 5,3 > secondSpatialDerivative_SourceHeadingDerivative;  ///< 5x3 matrix describing how the field spatial gradient changes with the source heading

    Eigen::Matrix<double, 9,3 > firstSpatialDerivative_secondSourceHeadingDerivative; ///< A list of 9x3 matricies describing how the field changes with the source heading [d(B*X)/dz; d(B*Y)/dz; d(B*Z)/dz]


    ScalorPotentialCalibrationJacobians();

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

/// @brief a class describing the field of a collection of scalor potentials
///
/// This class assumes the system can be described by spherical harmonic scalor
/// potentials at multiple locations.  See. "PUBLICATION HERE"
class ScalorPotential
{
public:
    struct srcCoeff {
        unsigned int order;
        double coeff;

        srcCoeff();
        srcCoeff(double value, unsigned int order);
    };

    /**
     * @brief The srcStruct struct describes the scalor potential for an individual source.
     *
     * The field given by a sources is the negative gradient of a scalor potential PHI according to the equation:
     * PHI = sum( (A_coeff(n)*r^n + B_coeff(n)*r^(-n-1))*P_n(cos(theta) , n= 1..inf)
     * where r is the distance from the point of interest to the source position. cos(theta) is the angle between the
     * r vector and the source direction. P_n is a legandre polynomial of order n. A_coeff and B_Coeff define the contributions
     * of each field shape.  Here we only keep the first terms in the summation.  Note: the length of A_Coeff and B_Coeff need not be the same.
     */
    class srcStruct {
    public:
        std::vector<srcCoeff> A_Coeff; /**< Ordered coefficients each associated with distances to an increasing positive power. **/
        std::vector<srcCoeff> B_Coeff; /**< Ordered coefficients each associated with distances to an increasing negative power. **/
        Eigen::Vector3d srcPosition;   /**< Source Position in meters **/
        Eigen::Vector3d srcDirection;  /**< Source Heading.  Should be a unit vector **/

        unsigned int getMaxOrder_A_Coeff() const;
        unsigned int getMaxOrder_B_Coeff() const;

        srcStruct();

        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    };

    /**
     * @brief ScalorPotential initializes an empty scalor potential, one with no sources.
     */
    ScalorPotential();

    ///
    /// \brief the constructor from a coil list.
    /// \param coilList The list of coils and their respective sources.
    /// \param dc_field_offset The dc offset field, if any.
    ///
    ScalorPotential( const std::vector<srcStruct>& srcList );

    ///
    /// \brief getState returns the state at the position requested
    /// \param position 3D position of interest
    /// \param sourceNumber The contribion limited to a specific source.  Default is -1, which is interperated as the total effect all sources
    /// \return The Scalor Potential State \see ScalorPotentialState
    ///
    ScalorPotentialState getState(const Eigen::Vector3d& position, int sourceNumber = -1 ) const;

    ///
    /// \brief getValue Returns the scalor potential magnitude at the position requested
    /// \param position 3D position of interest
    /// \param sourceNumber The contribution limited to a specific source.  Default is -1, which is interperated as the total effect of all sources
    /// \return The value of the scalor potential state
    ///
    double getValue( const Eigen::Vector3d& position, int sourceNumber = -1 ) const;

    ///
    /// \brief getGradient Returns the spacial gradient of the potential at the position requested
    /// \param position 3D position of interest
    /// \param sourceNumber The contribution limited to a specific source.  Default is -1, which is interperated as the total effect of all sources
    /// \return The value of the scalor potential spatial gradeint
    ///
    Eigen::Vector3d getGradient( const Eigen::Vector3d& position, int sourceNumber = -1 ) const;

    // These two functions convert a vector packing of the gradient to a matrix packing and vice versa
    ///@brief converts a 5x1 gradient vector into a symetric zero trace 3x3 matrix
    /// @param gradVector is the desired vector to be remaped
    static Eigen::Matrix3d remapSecondDerivativeVec(const Vector5d& gradVector );

    ///@brief converts a 3x3 gradient marix into the 5x1 gradient vector
    /// @param gradMatrix is the desired matrix to be remaped
    ///
    /// This function does not check to verify the matrix is indeed symetric and zero trace
    static Vector5d remapSecondDerivativeMat(const Eigen::Matrix3d& gradMatrix);


    ///@brief returns the number of sources for the given coil
    unsigned int getNumberOfSources( ) const;

    ///
    /// \brief getSourceStruct returns a copy ofthe source structure
    /// \param sourceNumber  The ordered number of the source to be copied.
    /// \return \see srcStruct.  If the number is out of bounds it will return an empty source.
    ///
    srcStruct getSourceStruct(unsigned int sourceNumber ) const;

    ///
    /// \brief setSourceStruct Sets a new source to replace an existing source.  Good to use with \see getSourceStruct for modifying a source configuration.
    /// \param sourceNumber The ordered number of the source to replace.
    /// \param newSrc The new source definition.
    ///
    /// This function will add empty sources until sourceNumber is reached if sourceNumber is greator than the current number of sources.
    ///
    void setSourceStruct(unsigned int sourceNumber, const srcStruct& newSrc);

    ///
    /// \brief removeSourceStruct  Removes a source from the list of sources.
    /// \param sourceNumber  The ordered number of the source to replace.
    ///
    void removeSourceStruct(unsigned int sourceNumber);

protected:
    std::vector<srcStruct> srcList;

    // functions for calculating field information
    static double LegandrePolynomial(double x, int order, int der = 0); ///< Calculates the Legandre Polynomial of a given order deritivative X <= (-1,1)
    static void srcFieldGradient(const Eigen::Vector3d& position, const srcStruct& src, ScalorPotentialState& currentState );

    // functions and data for helping with calibration
    ScalorPotentialCalibrationJacobians srcCalibrationInformation(const Eigen::Vector3d& position, unsigned int srcNum ) const; // generates the matricies required for calibration


    int numCalParameters;

    template <typename Derived>
    void packCalibrationState(Eigen::MatrixBase<Derived> const & stateVector_) const
    {
        assert( stateVector_.rows() == numCalParameters );

        /**< Returns the vector packing of the node state **/
        Eigen::MatrixBase<Derived>& stateVector = const_cast< Eigen::MatrixBase<Derived>& >(stateVector_);

        stateVector.setZero();
        int vectorIndex = 0;

        std::vector<srcStruct>::const_iterator srcIT = srcList.begin();

        for( ; srcIT != srcList.end(); srcIT ++ )
        {
            std::vector<srcCoeff>::const_iterator coeffIT = srcIT->A_Coeff.begin();

            // unpack A coeff
            for( ; coeffIT != srcIT->A_Coeff.end(); coeffIT ++, vectorIndex ++ )
            {
                stateVector(vectorIndex) = coeffIT->coeff;
            }

            // unpack B coeff
            for( coeffIT = srcIT->B_Coeff.begin(); coeffIT != srcIT->B_Coeff.end(); coeffIT ++, vectorIndex ++ )
            {
                stateVector(vectorIndex) = coeffIT->coeff;
            }

            // unpack direction
            for( int i=0; i<3; i++, vectorIndex ++ )
            {
                stateVector(vectorIndex) = srcIT->srcDirection(i);
            }

            // unpack position
            for( int i=0; i<3; i++, vectorIndex ++ )
            {
                stateVector(vectorIndex) = srcIT->srcPosition(i);
            }
        }

        return;
    }

    template <typename Derived>
    void unpackCalibrationState(Eigen::MatrixBase<Derived> const & stateVector_)
    {

        Eigen::MatrixBase<Derived>& stateVector = const_cast< Eigen::MatrixBase<Derived>& >(stateVector_);

        assert( ('State vector provided for unpacking is the worng size', stateVector.rows() == getNumCalibrationParameters() ));

        int vectorIndex = 0;

        std::vector<srcStruct>::iterator srcIT = srcList.begin();

        for( ; srcIT != srcList.end(); srcIT ++ )
        {
            std::vector<srcCoeff>::iterator coeffIT = srcIT->A_Coeff.begin();

            // unpack A coeff
            for( ; coeffIT != srcIT->A_Coeff.end(); coeffIT ++, vectorIndex ++ )
            {
                coeffIT->coeff = stateVector(vectorIndex);
            }

            // unpack B coeff
            for( coeffIT = srcIT->B_Coeff.begin(); coeffIT != srcIT->B_Coeff.end(); coeffIT ++, vectorIndex ++ )
            {
                coeffIT->coeff = stateVector(vectorIndex);
            }

            // unpack direction
            for( int i=0; i<3; i++, vectorIndex ++ )
            {
                srcIT->srcDirection(i) = stateVector(vectorIndex);
            }
            assert(('Direction has zero length!',srcIT->srcDirection.norm() != 0 ));
           // srcIT->srcDirection.normalize();

            // unpack position
            for( int i=0; i<3; i++, vectorIndex ++ )
            {
                srcIT->srcPosition(i) = stateVector(vectorIndex);
            }
        }

        assert( vectorIndex == getNumCalibrationParameters());
    }

    template <typename Derived, typename Derived2>
    void packJacobians(int srcNum, double current, const Eigen::Vector3d& pos, const Eigen::Vector3d& fieldError, Eigen::MatrixBase<Derived> const & firstTransposed_, Eigen::MatrixBase<Derived2> const & hessian_ )
    {

        Eigen::MatrixBase<Derived>& firstTransposed = const_cast< Eigen::MatrixBase<Derived>& >(firstTransposed_);
        Eigen::MatrixBase<Derived2>& hessian = const_cast< Eigen::MatrixBase<Derived2>& >(hessian_);


        int numParam = getNumCalibrationParameters(srcNum);
        assert( ( firstTransposed.rows() == numParam ));
        assert( (firstTransposed.cols() == 3 ));
        assert( hessian.rows() == numParam );
        assert( hessian.cols() == numParam );

        ScalorPotentialCalibrationJacobians jacobStruct =  srcCalibrationInformation(pos, srcNum );

        firstTransposed.setZero();

        int Ft_a_rowCount = 0;
        int numA = srcList[srcNum].A_Coeff.size();
        int numB = srcList[srcNum].B_Coeff.size();

        if( numA )
        {
            firstTransposed.block(Ft_a_rowCount,0,numA,3) = jacobStruct.firstSpatialDerivative_A_CoeffDerivative.transpose()*current;

            Eigen::Matrix3d dF_dAdH;
            for( int i=0; i< numA; i++ )
            {
                dF_dAdH.col(0) = jacobStruct.firstSpatialDerivative_dA_dHeading.block(0,i,3,1);
                dF_dAdH.col(1) = jacobStruct.firstSpatialDerivative_dA_dHeading.block(3,i,3,1);
                dF_dAdH.col(2) = jacobStruct.firstSpatialDerivative_dA_dHeading.block(6,i,3,1);

                hessian.block(Ft_a_rowCount+i,numA+numB,1,3) = fieldError.transpose()*dF_dAdH;
            }
            hessian.block(numA+numB, Ft_a_rowCount, 3,numA) = hessian.block(Ft_a_rowCount,numA+numB,numA,3).transpose();

//            hessian.block(Ft_a_rowCount,numA+numB,numA,3) = jacobStruct.firstSpatialDerivative_dA_dHeading.block(0,0,3,numA).transpose()*fieldError(0)
//                                                          + jacobStruct.firstSpatialDerivative_dA_dHeading.block(3,0,3,numA).transpose()*fieldError(1)
//                                                          + jacobStruct.firstSpatialDerivative_dA_dHeading.block(6,0,3,numA).transpose()*fieldError(2);
//            hessian.block(numA+numB, Ft_a_rowCount, 3,numA) = hessian.block(Ft_a_rowCount,numA+numB,numA,3).transpose();

            Eigen::Matrix3d dF_dAdP;
            for( int i=0; i< numA; i++ )
            {
                dF_dAdP.col(0) = jacobStruct.firstSpatialDerivative_dA_dPosition.block(0,i,3,1);
                dF_dAdP.col(1) = jacobStruct.firstSpatialDerivative_dA_dPosition.block(3,i,3,1);
                dF_dAdP.col(2) = jacobStruct.firstSpatialDerivative_dA_dPosition.block(6,i,3,1);

                hessian.block(Ft_a_rowCount+i,numA+numB+3,1,3) = fieldError.transpose()*dF_dAdP;
            }
            hessian.block(numA+numB+3, Ft_a_rowCount, 3,numA) = hessian.block(Ft_a_rowCount,numA+numB+3,numA,3).transpose();

//            hessian.block(Ft_a_rowCount,numA+numB+3,numA,3) = jacobStruct.firstSpatialDerivative_dA_dPosition.block(0,0,3,numA).transpose()*fieldError(0)
//                                                          +   jacobStruct.firstSpatialDerivative_dA_dPosition.block(3,0,3,numA).transpose()*fieldError(1)
//                                                          +   jacobStruct.firstSpatialDerivative_dA_dPosition.block(6,0,3,numA).transpose()*fieldError(2);
//            hessian.block(numA+numB+3, Ft_a_rowCount, 3,numA) = hessian.block(Ft_a_rowCount,numA+numB+3,numA,3).transpose();

            Ft_a_rowCount += numA;

        }
        if( numB )
        {
            firstTransposed.block(Ft_a_rowCount,0,numB,3) = jacobStruct.firstSpatialDerivative_B_CoeffDerivative.transpose()*current;

            Eigen::Matrix3d dF_dBdH;
            for( int i=0; i< numB; i++ )
            {
                dF_dBdH.col(0) = jacobStruct.firstSpatialDerivative_dB_dHeading.block(0,i,3,1);
                dF_dBdH.col(1) = jacobStruct.firstSpatialDerivative_dB_dHeading.block(3,i,3,1);
                dF_dBdH.col(2) = jacobStruct.firstSpatialDerivative_dB_dHeading.block(6,i,3,1);

                hessian.block(Ft_a_rowCount+i,numA+numB,1,3) = fieldError.transpose()*dF_dBdH;
            }
            hessian.block(numA+numB, Ft_a_rowCount, 3,numB) = hessian.block(Ft_a_rowCount,numA+numB,numB,3).transpose();

            Eigen::Matrix3d dF_dBdP;
            for( int i=0; i< numB; i++ )
            {
                dF_dBdP.col(0) = jacobStruct.firstSpatialDerivative_dB_dPosition.block(0,i,3,1);
                dF_dBdP.col(1) = jacobStruct.firstSpatialDerivative_dB_dPosition.block(3,i,3,1);
                dF_dBdP.col(2) = jacobStruct.firstSpatialDerivative_dB_dPosition.block(6,i,3,1);

                hessian.block(Ft_a_rowCount+i,numA+numB+3,1,3) = fieldError.transpose()*dF_dBdP;
            }
            hessian.block(numA+numB+3, Ft_a_rowCount, 3,numB) = hessian.block(Ft_a_rowCount,numA+numB+3,numB,3).transpose();


            Ft_a_rowCount += numB;
        }


        // heading terms
        firstTransposed.block(Ft_a_rowCount,0,3,3) = jacobStruct.firstSpatialDerivative_SourceHeadingDerivative.transpose()*current;

        //Eigen::Matrix3d dF_dHdH;
        for( int i=0; i< 3; i++ )
        {
            hessian.block(Ft_a_rowCount+i,numA+numB,1,3) += fieldError.transpose()*jacobStruct.firstSpatialDerivative_secondSourceHeadingDerivative.block(3*i,0,3,3);
        }

       /* hessian.block(Ft_a_rowCount,numA+numB,3,3) = jacobStruct.firstSpatialDerivative_secondSourceHeadingDerivative.block(0,0,3,3).transpose()*fieldError(0)
                                                   + jacobStruct.firstSpatialDerivative_secondSourceHeadingDerivative.block(3,0,3,3).transpose()*fieldError(1)
                                                   + jacobStruct.firstSpatialDerivative_secondSourceHeadingDerivative.block(6,0,3,3).transpose()*fieldError(2);

        hessian.block(numA+numB, Ft_a_rowCount, 3,3) = jacobStruct.firstSpatialDerivative_secondSourceHeadingDerivative.block(0,0,3,3)*fieldError(0)
                                                     + jacobStruct.firstSpatialDerivative_secondSourceHeadingDerivative.block(3,0,3,3)*fieldError(1)
                                                     + jacobStruct.firstSpatialDerivative_secondSourceHeadingDerivative.block(6,0,3,3)*fieldError(2);
       */
        for( int i=0; i< 3; i++ )
        {
            hessian.block(Ft_a_rowCount+i,numA+numB+3,1,3) += fieldError.transpose()*remapSecondDerivativeVec(jacobStruct.secondSpatialDerivative_SourceHeadingDerivative.col(i));
        }
        hessian.block(numA+numB+3,Ft_a_rowCount,3,3) = hessian.block(Ft_a_rowCount,numA+numB+3,3,3).transpose();

//        hessian.block(Ft_a_rowCount,numA+numB+3,3,3) = remapSecondDerivativeVec(jacobStruct.secondSpatialDerivative_SourceHeadingDerivative.col(0))*fieldError(0)
//                                                   +   remapSecondDerivativeVec(jacobStruct.secondSpatialDerivative_SourceHeadingDerivative.col(1))*fieldError(1)
//                                                   +   remapSecondDerivativeVec(jacobStruct.secondSpatialDerivative_SourceHeadingDerivative.col(2))*fieldError(2);
//        hessian.block(numA+numB+3, Ft_a_rowCount, 3,3) = (remapSecondDerivativeVec(jacobStruct.secondSpatialDerivative_SourceHeadingDerivative.col(0))*fieldError(0)
//                                                      +   remapSecondDerivativeVec(jacobStruct.secondSpatialDerivative_SourceHeadingDerivative.col(1))*fieldError(1)
//                                                      +   remapSecondDerivativeVec(jacobStruct.secondSpatialDerivative_SourceHeadingDerivative.col(2))*fieldError(2)).transpose();


        Ft_a_rowCount += 3;

        // position terms
        firstTransposed.block(Ft_a_rowCount,0,3,3) = jacobStruct.firstSpatialDerivative_SourcePositionDerivative*current;

        hessian.block(Ft_a_rowCount,numA+numB+3,3,3) = remapSecondDerivativeVec(jacobStruct.secondSpatialDerivative_SourcePositionDerivative.col(0))*fieldError(0)
                                                   +   remapSecondDerivativeVec(jacobStruct.secondSpatialDerivative_SourcePositionDerivative.col(1))*fieldError(1)
                                                   +   remapSecondDerivativeVec(jacobStruct.secondSpatialDerivative_SourcePositionDerivative.col(2))*fieldError(2);

        Ft_a_rowCount += 3;

        hessian *= current;


    }

    int getNumCalibrationParameters(int srcNum = -1) const;

    friend class ElectromagnetCalibration;

};

#endif // ScalorPotential_H
